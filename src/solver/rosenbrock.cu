// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <chrono>
#include <iostream>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/util/cuda_param.hpp>
#include <vector>

#include "cublas_v2.h"

namespace micm
{
  namespace cuda
  {
    /// CUDA kernel to compute alpha - J[i] for each element i at the diagnoal of matrix J
    __global__ void AlphaMinusJacobianKernel(double* d_jacobian, const double alpha, CudaRosenbrockSolverParam devstruct)
    {
      // Global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local variables
      size_t quotient, index_as_remainder;
      const size_t num_diagonal_elements = devstruct.jacobian_diagonal_elements_size_;
      const size_t num_grid_cells = devstruct.num_grid_cells_;

      if (tid < num_grid_cells * num_diagonal_elements)
      {
        quotient = tid / num_diagonal_elements;
        index_as_remainder = tid - num_diagonal_elements * quotient;  // % operator may be more expensive
        d_jacobian[devstruct.jacobian_diagonal_elements_[index_as_remainder] + quotient] += alpha;
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

      /// Create a struct whose members contain the addresses in the device memory.
      CudaRosenbrockSolverParam devstruct;
      cudaMalloc(&(devstruct.errors_input_), errors_bytes);
      cudaMalloc(&(devstruct.errors_output_), errors_bytes);
      cudaMalloc(&(devstruct.jacobian_diagonal_elements_), jacobian_diagonal_elements_bytes);

      /// Copy the data from host to device
      cudaMemcpy(
          devstruct.jacobian_diagonal_elements_,
          hoststruct.jacobian_diagonal_elements_,
          jacobian_diagonal_elements_bytes,
          cudaMemcpyHostToDevice);

      devstruct.num_grid_cells_ = hoststruct.num_grid_cells_;
      devstruct.errors_size_ = hoststruct.errors_size_;
      devstruct.jacobian_diagonal_elements_size_ = hoststruct.jacobian_diagonal_elements_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members and temporary variables of class "CudaLuDecomposition" on the device
    void FreeConstData(CudaRosenbrockSolverParam& devstruct)
    {
      cudaFree(devstruct.errors_input_);
      cudaFree(devstruct.errors_output_);
      cudaFree(devstruct.jacobian_diagonal_elements_);
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
        const CudaVectorMatrixParam y_old_param,
        const CudaVectorMatrixParam y_new_param,
        const RosenbrockSolverParameters ros_param,
        CudaRosenbrockSolverParam devstruct,
        const size_t n,
        bool is_first_call)
    {
      double* d_y_old = y_old_param.d_data_;
      double* d_y_new = y_new_param.d_data_;
      double* d_errors_input = devstruct.errors_input_;
      double* d_errors_output = devstruct.errors_output_;
      const double atol = ros_param.absolute_tolerance_;
      const double rtol = ros_param.relative_tolerance_;

      // Declares a dynamically-sized shared memory array.
      // The size of this array is determined at runtime when the kernel is launched.
      // Shared memory is shared among all threads within the same block.
      extern __shared__ double sdata[];

      // Local thread ID within a threadblock
      size_t l_tid = threadIdx.x;

      // Global thread ID
      size_t g_tid = blockIdx.x * (BLOCK_SIZE * 2) + threadIdx.x;

      if (is_first_call)
      {
        // Temporary device variables
        double d_ymax, d_scale;

        // Load two elements by one thread and do first add of reduction
        sdata[l_tid] = 0.0;
        for (int i = 0; i < 2; ++i)
        {
          if (g_tid < n)
          {
            d_ymax = max(fabs(d_y_old[g_tid]), fabs(d_y_new[g_tid]));
            d_scale = atol + rtol * d_ymax;
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
        const CudaVectorMatrixParam y_old_param,
        const CudaVectorMatrixParam y_new_param,
        const RosenbrockSolverParameters ros_param,
        CudaRosenbrockSolverParam devstruct)
    {
      // Temporary device variables
      double d_ymax, d_scale;
      double* d_y_old = y_old_param.d_data_;
      double* d_y_new = y_new_param.d_data_;
      double* d_errors = devstruct.errors_input_;
      double atol = ros_param.absolute_tolerance_;
      double rtol = ros_param.relative_tolerance_;
      const size_t num_elements = devstruct.errors_size_;

      // Global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid < num_elements)
      {
        d_ymax = max(fabs(d_y_old[tid]), fabs(d_y_new[tid]));
        d_scale = atol + rtol * d_ymax;
        d_errors[tid] = d_errors[tid] / d_scale;
      }
    }

    // Host code that will launch the AlphaMinusJacobian CUDA kernel
    void AlphaMinusJacobianDriver(
        double* h_jacobian,
        const size_t num_elements,
        const double alpha,
        const CudaRosenbrockSolverParam& devstruct)
    {
      // device pointers (will not be needed after adding the CudaSparseMatrix class)
      double* d_jacobian;
      cudaMalloc(&d_jacobian, sizeof(double) * num_elements);
      cudaMemcpy(d_jacobian, h_jacobian, sizeof(double) * num_elements, cudaMemcpyHostToDevice);

      // kernel call
      size_t num_blocks =
          (devstruct.jacobian_diagonal_elements_size_ * devstruct.num_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      AlphaMinusJacobianKernel<<<num_blocks, BLOCK_SIZE>>>(d_jacobian, alpha, devstruct);

      cudaDeviceSynchronize();
      cudaMemcpy(h_jacobian, d_jacobian, sizeof(double) * num_elements, cudaMemcpyDeviceToHost);
      cudaFree(d_jacobian);
    }

    // Host code that will launch the NormalizedError CUDA kernel
    double NormalizedErrorDriver(
        const CudaVectorMatrixParam& y_old_param,
        const CudaVectorMatrixParam& y_new_param,
        const CudaVectorMatrixParam& errors_param,
        const RosenbrockSolverParameters& ros_param,
        cublasHandle_t handle,
        CudaRosenbrockSolverParam devstruct)
    {
      double normalized_error;
      const size_t num_elements = devstruct.errors_size_;

      if (devstruct.errors_size_ != errors_param.number_of_elements_)
      {
        throw std::runtime_error("devstruct.errors_input_ and errors_param have different sizes.");
      }
      cudaError_t err =
          cudaMemcpy(devstruct.errors_input_, errors_param.d_data_, sizeof(double) * num_elements, cudaMemcpyDeviceToDevice);

      if (num_elements > 1000000)
      {
        // call cublas APIs
        size_t num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;
        ScaledErrorKernel<<<num_blocks, BLOCK_SIZE>>>(y_old_param, y_new_param, ros_param, devstruct);
        // call cublas function to perform the norm:
        // https://docs.nvidia.com/cuda/cublas/index.html?highlight=dnrm2#cublas-t-nrm2
        cublasStatus_t stat = cublasDnrm2(handle, num_elements, devstruct.errors_input_, 1, &normalized_error);
        if (stat != CUBLAS_STATUS_SUCCESS)
        {
          printf(cublasGetStatusString(stat));
          throw std::runtime_error("Error while calling cublasDnrm2.");
        }
        normalized_error = normalized_error * std::sqrt(1.0 / num_elements);
      }
      else
      {
        // call CUDA implementation
        size_t num_blocks = std::ceil(std::ceil(num_elements * 1.0 / BLOCK_SIZE) / 2.0);
        num_blocks = num_blocks < 1 ? 1 : num_blocks;
        size_t new_blocks;
        bool is_first_call;

        is_first_call = true;
        // Kernel call
        NormalizedErrorKernel<<<num_blocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
            y_old_param, y_new_param, ros_param, devstruct, num_elements, is_first_call);
        is_first_call = false;
        while (num_blocks > 1)
        {
          std::swap(devstruct.errors_input_, devstruct.errors_output_);
          // Update grid size
          new_blocks = std::ceil(std::ceil(num_blocks * 1.0 / BLOCK_SIZE) / 2.0);
          if (new_blocks <= 1)
          {
            NormalizedErrorKernel<<<1, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
                y_old_param, y_new_param, ros_param, devstruct, num_blocks, is_first_call);
            break;
          }
          NormalizedErrorKernel<<<new_blocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double)>>>(
              y_old_param, y_new_param, ros_param, devstruct, num_blocks, is_first_call);
          num_blocks = new_blocks;
        }
        cudaDeviceSynchronize();

        cudaMemcpy(&normalized_error, &devstruct.errors_output_[0], sizeof(double), cudaMemcpyDeviceToHost);
        normalized_error = std::sqrt(normalized_error / num_elements);
      }  // end of if-else for CUDA/CUBLAS implementation
      return std::max(normalized_error, 1.0e-10);
    }  // end of NormalizedErrorDriver function
  }    // namespace cuda
}  // namespace micm
