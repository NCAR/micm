// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/types.hpp>

#include <cublas_v2.h>

#include <type_traits>

namespace micm::cuda
{
  /// CUDA kernel to compute alpha - J[i] for each element i at the diagonal of Jacobian matrix
  __global__ void AlphaMinusJacobianKernel(
      CudaMatrixParam jacobian_param,
      const Real alpha,
      const CudaJacobianDiagonalElementsParam jacobian_diagonal_elements_param)
  {
    // Calculate global thread ID
    const Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    // Local device variables
    Real* d_jacobian = jacobian_param.d_data_;
    Index quotient, index_as_remainder;
    const Index number_of_diagonal_elements = jacobian_diagonal_elements_param.size_;
    const Index number_of_grid_cells = jacobian_param.number_of_grid_cells_;
    const Index cuda_matrix_vector_length = jacobian_param.vector_length_;
    const Index number_of_groups = (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
    const Index number_of_diagonal_elements_per_group = number_of_diagonal_elements * cuda_matrix_vector_length;
    const Index number_of_non_zeros_per_group =
        jacobian_param.number_of_elements_ / (number_of_groups * cuda_matrix_vector_length);
    const Index total_number_of_diagonal_elements = number_of_diagonal_elements_per_group * number_of_groups;
    const Index group_id = tid / number_of_diagonal_elements_per_group;
    const Index local_tid = tid - group_id * number_of_diagonal_elements_per_group;

    // Shift the index for different groups
    d_jacobian += group_id * number_of_non_zeros_per_group * cuda_matrix_vector_length;

    if (tid < total_number_of_diagonal_elements)
    {
      quotient = local_tid / cuda_matrix_vector_length;
      index_as_remainder = local_tid - cuda_matrix_vector_length * quotient;  // % operator may be more expensive
      d_jacobian[jacobian_diagonal_elements_param.data_[quotient] + index_as_remainder] += alpha;
    }
  }

  // Specific CUDA device function to do reduction within a warp
  // Use volatile to prevent compiler optimization (caching in registers)
  // No need to synchronize threads in the same warp
  __device__ void WarpReduce(volatile Real* sdata, const Index tid)
  {
    if (BLOCK_SIZE >= 64)
    {
      sdata[tid] = sdata[tid] + sdata[tid + 32];
    }
    sdata[tid] = sdata[tid] + sdata[tid + 16];
    sdata[tid] = sdata[tid] + sdata[tid + 8];
    sdata[tid] = sdata[tid] + sdata[tid + 4];
    sdata[tid] = sdata[tid] + sdata[tid + 2];
    sdata[tid] = sdata[tid] + sdata[tid + 1];
  }

  // CUDA kernel to compute the scaled norm of the vector errors; CUDA kernel does not take reference as argument
  // Modified version from NVIDIA's reduction example:
  // https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
  __global__ void NormalizedErrorKernel(
      const CudaMatrixParam y_old_param,
      const CudaMatrixParam y_new_param,
      const CudaMatrixParam absolute_tolerance_param,
      const Real relative_tolerance,
      const Index number_of_elements,
      bool is_first_call,
      const CudaErrorParam errors_calc_param)
  {
    Real* const d_errors_input = errors_calc_param.errors_input_;
    Real* const d_errors_output = errors_calc_param.errors_output_;

    // Declares a dynamically-sized shared memory array.
    // The size of this array is determined at runtime when the kernel is launched.
    // Shared memory is shared among all threads within the same block.
    extern __shared__ Real sdata[];

    // Calculate local thread ID within a threadblock
    const Index l_tid = threadIdx.x;

    // Calculate global thread ID
    Index g_tid = blockIdx.x * (BLOCK_SIZE * 2) + threadIdx.x;

    if (is_first_call)
    {
      // Load two elements by one thread and do first add of reduction
      sdata[l_tid] = 0.0;
      for (Index i = 0; i < 2; ++i)
      {
        if (g_tid < number_of_elements)
        {
          d_errors_input[g_tid] = d_errors_input[g_tid] * d_errors_input[g_tid];
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
      if (g_tid < number_of_elements)
      {
        sdata[l_tid] += d_errors_input[g_tid];
      }
      g_tid += BLOCK_SIZE;
      if (g_tid < number_of_elements)
      {
        sdata[l_tid] += d_errors_input[g_tid];
      }
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
    {
      WarpReduce(sdata, l_tid);
    }

    // Let the thread 0 of this threadblock write its result to output array, inexed by this threadblock
    if (l_tid == 0)
    {
      d_errors_output[blockIdx.x] = sdata[0];
    }
  }

  // CUDA kernel to compute the scaled vectors
  __global__ void ScaledErrorKernel(
      const CudaMatrixParam y_old_param,
      const CudaMatrixParam y_new_param,
      const CudaMatrixParam absolute_tolerance_param,
      const Real relative_tolerance,
      const CudaErrorParam errors_param)
  {
    // Local device variables
    Real d_ymax, d_scale;
    const Real* const d_y_old = y_old_param.d_data_;
    const Real* const d_y_new = y_new_param.d_data_;
    Real* const d_errors = errors_param.errors_input_;
    const Real* const atol = absolute_tolerance_param.d_data_;
    const Index number_of_elements = errors_param.errors_size_;
    const Index number_of_grid_cells = y_old_param.number_of_grid_cells_;
    const Index cuda_matrix_vector_length = y_old_param.vector_length_;
    const Index number_of_variables = absolute_tolerance_param.number_of_elements_;
    const Index full_group =
        (number_of_grid_cells / cuda_matrix_vector_length) * cuda_matrix_vector_length * number_of_variables;

    // Calculate global thread ID
    const Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

    // Calculate the index for absolute tolerance array
    const Index atol_idx = (tid / cuda_matrix_vector_length) % number_of_variables;

    // compute the elements over the groups which fit exactly into the vector length parameter
    if (tid < full_group)
    {
      d_ymax = max(fabs(d_y_old[tid]), fabs(d_y_new[tid]));
      d_scale = atol[atol_idx] + relative_tolerance * d_ymax;
      d_errors[tid] = d_errors[tid] / d_scale;
    }

    // compute the elements over the remaining group that may be partially filled
    const Index partial_group = number_of_grid_cells % cuda_matrix_vector_length;
    if (partial_group > 0 && tid >= full_group && tid < number_of_elements)
    {
      const Index row_idx = tid % cuda_matrix_vector_length;
      const Index local_tid = full_group + atol_idx * cuda_matrix_vector_length + row_idx;
      if (row_idx < partial_group)
      {
        d_ymax = max(fabs(d_y_old[local_tid]), fabs(d_y_new[local_tid]));
        d_scale = atol[atol_idx] + relative_tolerance * d_ymax;
        d_errors[local_tid] = d_errors[local_tid] / d_scale;
      }
      else
      {
        d_errors[local_tid] = 0.0;
      }
    }
  }

  // Host code that will launch the AlphaMinusJacobian CUDA kernel
  void AlphaMinusJacobianDriver(
      CudaMatrixParam& jacobian_param,
      const Real& alpha,
      const CudaJacobianDiagonalElementsParam& jacobian_diagonal_elements_param)
  {
    // We will add alpha to the padding elements as well to simplify the kernel code
    const Index number_of_groups =
        (jacobian_param.number_of_grid_cells_ + jacobian_param.vector_length_ - 1) / jacobian_param.vector_length_;
    const Index number_of_blocks =
        (jacobian_diagonal_elements_param.size_ * number_of_groups * jacobian_param.vector_length_ + BLOCK_SIZE - 1) /
        BLOCK_SIZE;
    AlphaMinusJacobianKernel<<<
        number_of_blocks,
        BLOCK_SIZE,
        0,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
        jacobian_param, alpha, jacobian_diagonal_elements_param);
  }

  // Host code that will launch the NormalizedError CUDA kernel
  Real NormalizedErrorDriver(
      const CudaMatrixParam& y_old_param,
      const CudaMatrixParam& y_new_param,
      const CudaMatrixParam& y_error_param,
      const CudaMatrixParam& absolute_tolerance_param,
      const Real relative_tolerance,
      CudaErrorParam& errors_param)
  {
    Real normalized_error;
    const Index number_of_grid_cells = y_old_param.number_of_grid_cells_;
    const Index number_of_species = absolute_tolerance_param.number_of_elements_;
    const Index number_of_elements = errors_param.errors_size_;

    if (number_of_elements != y_error_param.number_of_elements_)
    {
      std::string msg = "mismatch in normalized error arrays. Expected: " + std::to_string(number_of_elements) +
                        " but got: " + std::to_string(y_error_param.number_of_elements_);
      throw micm::MicmException(MICM_ERROR_CATEGORY_INTERNAL, MICM_INTERNAL_ERROR_CODE_GENERAL, msg);
    }
    CHECK_CUDA_ERROR(
        cudaMemcpyAsync(
            errors_param.errors_input_,
            y_error_param.d_data_,
            sizeof(Real) * number_of_elements,
            cudaMemcpyDeviceToDevice,
            micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
        "cudaMemcpy");

    // Launch the ScaledErrorKernel to compute the scaled errors
    Index number_of_blocks = (number_of_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;
    ScaledErrorKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
        y_old_param, y_new_param, absolute_tolerance_param, relative_tolerance, errors_param);

    if (number_of_elements > 1000000)
    {
      // call cublas function to perform the norm:
      // https://docs.nvidia.com/cuda/cublas/index.html?highlight=dnrm2#cublas-t-nrm2
      static_assert(std::is_same_v<micm::Real, double>, "cuBLAS D-routines require Real == double");
      CHECK_CUBLAS_ERROR(
          cublasDnrm2(micm::cuda::GetCublasHandle(), number_of_elements, errors_param.errors_input_, 1, &normalized_error),
          "cublasDnrm2");
      cudaStreamSynchronize(micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
      normalized_error = normalized_error * std::sqrt(1.0 / (number_of_grid_cells * number_of_species));
    }
    else
    {
      // call CUDA implementation
      number_of_blocks = (number_of_elements + 2 * BLOCK_SIZE - 1) / (2 * BLOCK_SIZE);
      number_of_blocks = number_of_blocks < 1 ? 1 : number_of_blocks;
      Index new_number_of_blocks;
      bool is_first_call = true;

      // Kernel call
      NormalizedErrorKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          BLOCK_SIZE * sizeof(Real),
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          y_old_param,
          y_new_param,
          absolute_tolerance_param,
          relative_tolerance,
          number_of_elements,
          is_first_call,
          errors_param);
      is_first_call = false;
      while (number_of_blocks > 1)
      {
        std::swap(errors_param.errors_input_, errors_param.errors_output_);
        // Update grid size
        new_number_of_blocks = (number_of_blocks + 2 * BLOCK_SIZE - 1) / (2 * BLOCK_SIZE);
        if (new_number_of_blocks <= 1)
        {
          NormalizedErrorKernel<<<
              1,
              BLOCK_SIZE,
              BLOCK_SIZE * sizeof(Real),
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
              y_old_param,
              y_new_param,
              absolute_tolerance_param,
              relative_tolerance,
              number_of_blocks,
              is_first_call,
              errors_param);
          break;
        }
        NormalizedErrorKernel<<<
            new_number_of_blocks,
            BLOCK_SIZE,
            BLOCK_SIZE * sizeof(Real),
            micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
            y_old_param,
            y_new_param,
            absolute_tolerance_param,
            relative_tolerance,
            number_of_blocks,
            is_first_call,
            errors_param);
        number_of_blocks = new_number_of_blocks;
      }

      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              &normalized_error,
              &errors_param.errors_output_[0],
              sizeof(Real),
              cudaMemcpyDeviceToHost,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      cudaStreamSynchronize(micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
      normalized_error = std::sqrt(normalized_error / (number_of_grid_cells * number_of_species));
    }  // end of if-else for CUDA/CUBLAS implementation
    const Real error_min = 1.0e-10;
    return std::max(normalized_error, error_min);
  }  // end of NormalizedErrorDriver function
}  // namespace micm::cuda
