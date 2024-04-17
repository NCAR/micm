// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <iostream>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/util/cuda_param.hpp>

#include "cublas_v2.h"

namespace micm
{
  namespace cuda
  {
    /// CUDA kernel to compute alpha - J[i] for each element i at the diagnoal of Jacobian matrix
    __global__ void
    AlphaMinusJacobianKernel(CudaMatrixParam jacobian_param, const double alpha, const CudaRosenbrockSolverParam devstruct)
    {
      // Calculate global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      double* d_jacobian = jacobian_param.d_data_;
      size_t quotient, index_as_remainder;
      const size_t number_of_diagonal_elements = devstruct.jacobian_diagonal_elements_size_;
      const size_t number_of_grid_cells = jacobian_param.number_of_grid_cells_;

      if (tid < number_of_grid_cells * number_of_diagonal_elements)
      {
        quotient = tid / number_of_grid_cells;
        index_as_remainder = tid - number_of_grid_cells * quotient;  // % operator may be more expensive
        d_jacobian[devstruct.jacobian_diagonal_elements_[quotient] + index_as_remainder] += alpha;
      }
    }

    /// This is the function that will allocate device memory
    ///   and copy const data for data members of class "CudaRosenbrockSolverParam"
    CudaRosenbrockSolverParam CopyConstData(CudaRosenbrockSolverParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t jacobian_diagonal_elements_bytes = sizeof(size_t) * hoststruct.jacobian_diagonal_elements_size_;

      /// Calculate the memory space of each temporary variable
      size_t errors_bytes = sizeof(double) * hoststruct.errors_size_;
      size_t tolerance_bytes = sizeof(double) * hoststruct.absolute_tolerance_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      CudaRosenbrockSolverParam devstruct;
      cudaMalloc(&(devstruct.errors_input_), errors_bytes);
      cudaMalloc(&(devstruct.errors_output_), errors_bytes);
      cudaMalloc(&(devstruct.jacobian_diagonal_elements_), jacobian_diagonal_elements_bytes);
      cudaMalloc(&(devstruct.absolute_tolerance_), tolerance_bytes);

      /// Copy the data from host to device
      cudaMemcpy(
          devstruct.jacobian_diagonal_elements_,
          hoststruct.jacobian_diagonal_elements_,
          jacobian_diagonal_elements_bytes,
          cudaMemcpyHostToDevice);
      
      cudaMemcpy(
          devstruct.absolute_tolerance_,
          hoststruct.absolute_tolerance_,
          tolerance_bytes,
          cudaMemcpyHostToDevice);

      devstruct.errors_size_ = hoststruct.errors_size_;
      devstruct.jacobian_diagonal_elements_size_ = hoststruct.jacobian_diagonal_elements_size_;
      devstruct.absolute_tolerance_size_ = hoststruct.absolute_tolerance_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members and temporary variables of class "CudaLuDecomposition" on the device
    void FreeConstData(CudaRosenbrockSolverParam& devstruct)
    {
      cudaFree(devstruct.errors_input_);
      cudaFree(devstruct.errors_output_);
      cudaFree(devstruct.jacobian_diagonal_elements_);
      cudaFree(devstruct.absolute_tolerance_);
    }

    // Specific CUDA device function to do reduction within a warp
    // Use volatile to prevent compiler optimization (caching in registers)
    // No need to synchronize threads in the same warp
    __device__ void warpReduce(volatile double* sdata, size_t tid)
    {
      if (BLOCK_SIZE >= 64)
        sdata[tid] += sdata[tid + 32];
      sdata[tid] += sdata[tid + 16];
      sdata[tid] += sdata[tid + 8];
      sdata[tid] += sdata[tid + 4];
      sdata[tid] += sdata[tid + 2];
      sdata[tid] += sdata[tid + 1];
    }

    // CUDA kernel to compute the scaled norm of the vector errors; CUDA kernel does not take reference as argument
    // Modified version from NVIDIA's reduction example:
    // https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
    __global__ void NormalizedErrorKernel(
        const CudaMatrixParam y_old_param,
        const CudaMatrixParam y_new_param,
        const RosenbrockSolverParameters ros_param,
        CudaRosenbrockSolverParam devstruct,
        const size_t n,
        bool is_first_call)
    {
      double* d_y_old = y_old_param.d_data_;
      double* d_y_new = y_new_param.d_data_;
      double* d_errors_input = devstruct.errors_input_;
      double* d_errors_output = devstruct.errors_output_;
      const double* atol = devstruct.absolute_tolerance_;
      const double rtol = ros_param.relative_tolerance_;
      const size_t number_of_species = y_old_param.number_of_elements_;

      // Declares a dynamically-sized shared memory array.
      // The size of this array is determined at runtime when the kernel is launched.
      // Shared memory is shared among all threads within the same block.
      extern __shared__ double sdata[];

      // Calculate local thread ID within a threadblock
      size_t l_tid = threadIdx.x;

      // Calculate global thread ID
      size_t g_tid = blockIdx.x * (BLOCK_SIZE * 2) + threadIdx.x;

      if (is_first_call)
      {
        // Local device variables
        double d_ymax, d_scale;

        // Load two elements by one thread and do first add of reduction
        sdata[l_tid] = 0.0;
        for (int i = 0; i < 2; ++i)
        {
          if (g_tid < n)
          {
            d_ymax = max(fabs(d_y_old[g_tid]), fabs(d_y_new[g_tid]));
            d_scale = atol[g_tid % number_of_species] + rtol * d_ymax;
            d_errors_input[g_tid] = d_errors_input[g_tid] * d_errors_input[g_tid] / (d_scale * d_scale);
            sdata[l_tid] += d_errors_input[g_tid];
          }
          g_tid += BLOCK_SIZE;
        }
        __syncthreads();
      }
      else
      {
        // Load two elements by one thread and do first add of reduction
        // Access the d_errors array directly if it is not the first call
        sdata[l_tid] = 0.0;
        if (g_tid < n)
          sdata[l_tid] += d_errors_input[g_tid];
        g_tid += BLOCK_SIZE;
        if (g_tid < n)
          sdata[l_tid] += d_errors_input[g_tid];
        __syncthreads();
      }

      // Start at 1/2 block stride, do the add, and divide by two each iteration
      if (BLOCK_SIZE >= 1024)
      {
        if (l_tid < 512)
        {
          sdata[l_tid] += sdata[l_tid + 512];
        }
        __syncthreads();
      }
      if (BLOCK_SIZE >= 512)
      {
        if (l_tid < 256)
        {
          sdata[l_tid] += sdata[l_tid + 256];
        }
        __syncthreads();
      }
      if (BLOCK_SIZE >= 256)
      {
        if (l_tid < 128)
        {
          sdata[l_tid] += sdata[l_tid + 128];
        }
        __syncthreads();
      }
      if (BLOCK_SIZE >= 128)
      {
        if (l_tid < 64)
        {
          sdata[l_tid] += sdata[l_tid + 64];
        }
        __syncthreads();
      }
      if (l_tid < 32)
        warpReduce(sdata, l_tid);

      // Let the thread 0 of this threadblock write its result to output array, inexed by this threadblock
      if (l_tid == 0)
        d_errors_output[blockIdx.x] = sdata[0];
    }

    // CUDA kernel to compute the scaled vectors; prepare the input for cublas call later
    __global__ void ScaledErrorKernel(
        const CudaMatrixParam y_old_param,
        const CudaMatrixParam y_new_param,
        const RosenbrockSolverParameters ros_param,
        CudaRosenbrockSolverParam devstruct)
    {
      // Local device variables
      double d_ymax, d_scale;
      double* d_y_old = y_old_param.d_data_;
      double* d_y_new = y_new_param.d_data_;
      double* d_errors = devstruct.errors_input_;
      double* atol = devstruct.absolute_tolerance_;
      double rtol = ros_param.relative_tolerance_;
      const size_t num_elements = devstruct.errors_size_;
      const size_t number_of_species = y_old_param.number_of_elements_;

      // Calculate global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid < num_elements)
      {
        d_ymax = max(fabs(d_y_old[tid]), fabs(d_y_new[tid]));
        d_scale = atol[tid % number_of_species] + rtol * d_ymax;
        d_errors[tid] = d_errors[tid] / d_scale;
      }
    }

    // Host code that will launch the AlphaMinusJacobian CUDA kernel
    void AlphaMinusJacobianDriver(
        CudaMatrixParam& jacobian_param,
        const double& alpha,
        const CudaRosenbrockSolverParam& devstruct)
    {
      size_t number_of_blocks =
          (devstruct.jacobian_diagonal_elements_size_ * jacobian_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      AlphaMinusJacobianKernel<<<number_of_blocks, BLOCK_SIZE>>>(jacobian_param, alpha, devstruct);
      cudaDeviceSynchronize();
    }

    // Host code that will launch the NormalizedError CUDA kernel
    double NormalizedErrorDriver(
        const CudaMatrixParam& y_old_param,
        const CudaMatrixParam& y_new_param,
        const CudaMatrixParam& errors_param,
        const RosenbrockSolverParameters& ros_param,
        cublasHandle_t handle,
        CudaRosenbrockSolverParam devstruct)
    {
      double normalized_error;
      const size_t number_of_elements = devstruct.errors_size_;

      if (number_of_elements != errors_param.number_of_elements_)
      {
        throw std::runtime_error("devstruct.errors_input_ and errors_param have different sizes.");
      }
      cudaError_t err = cudaMemcpy(
          devstruct.errors_input_, errors_param.d_data_, sizeof(double) * number_of_elements, cudaMemcpyDeviceToDevice);

      if (number_of_elements > 1000000)
      {
        // call cublas APIs
        size_t number_of_blocks = (number_of_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;
        ScaledErrorKernel<<<number_of_blocks, BLOCK_SIZE>>>(y_old_param, y_new_param, ros_param, devstruct);
        // call cublas function to perform the norm:
        // https://docs.nvidia.com/cuda/cublas/index.html?highlight=dnrm2#cublas-t-nrm2
        cublasStatus_t stat = cublasDnrm2(handle, number_of_elements, devstruct.errors_input_, 1, &normalized_error);
        if (stat != CUBLAS_STATUS_SUCCESS)
        {
          printf(cublasGetStatusString(stat));
          throw std::runtime_error("Error while calling cublasDnrm2.");
        }
        normalized_error = normalized_error * std::sqrt(1.0 / number_of_elements);
      }
      else
      {
        // call CUDA implementation
        size_t number_of_blocks = std::ceil(std::ceil(number_of_elements * 1.0 / BLOCK_SIZE) / 2.0);
        number_of_blocks = number_of_blocks < 1 ? 1 : number_of_blocks;
        size_t new_number_of_blocks;
        bool is_first_call = true;

        // Kernel call
        NormalizedErrorKernel<<<number_of_blocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
            y_old_param, y_new_param, ros_param, devstruct, number_of_elements, is_first_call);
        is_first_call = false;
        while (number_of_blocks > 1)
        {
          std::swap(devstruct.errors_input_, devstruct.errors_output_);
          // Update grid size
          new_number_of_blocks = std::ceil(std::ceil(number_of_blocks * 1.0 / BLOCK_SIZE) / 2.0);
          if (new_number_of_blocks <= 1)
          {
            NormalizedErrorKernel<<<1, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
                y_old_param, y_new_param, ros_param, devstruct, number_of_blocks, is_first_call);
            break;
          }
          NormalizedErrorKernel<<<new_number_of_blocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
              y_old_param, y_new_param, ros_param, devstruct, number_of_blocks, is_first_call);
          number_of_blocks = new_number_of_blocks;
        }
        cudaDeviceSynchronize();

        cudaMemcpy(&normalized_error, &devstruct.errors_output_[0], sizeof(double), cudaMemcpyDeviceToHost);
        normalized_error = std::sqrt(normalized_error / number_of_elements);
      }  // end of if-else for CUDA/CUBLAS implementation
      return std::max(normalized_error, 1.0e-10);
    }  // end of NormalizedErrorDriver function
  }    // namespace cuda
}  // namespace micm
