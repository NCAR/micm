// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <chrono>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs the "solve" function on the device
    __global__ void
    SolveKernel(CudaMatrixParam x_param, const CudaMatrixParam ALU_param, const LinearSolverInPlaceParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::size_t* const __restrict__ d_nLij = devstruct.nLij_;
      const std::pair<std::size_t, std::size_t>* __restrict__ d_Lij_yj = devstruct.Lij_yj_;
      const std::pair<std::size_t, std::size_t>* const __restrict__ d_nUij_Uii = devstruct.nUij_Uii_;
      const std::pair<std::size_t, std::size_t>* __restrict__ d_Uij_xj = devstruct.Uij_xj_;
      const std::size_t d_nLij_size = devstruct.nLij_size_;
      const std::size_t d_nUij_Uii_size = devstruct.nUij_Uii_size_;

      double* __restrict__ d_ALU = ALU_param.d_data_;
      double* d_x = x_param.d_data_;
      const std::size_t number_of_grid_cells = x_param.number_of_grid_cells_;
      const std::size_t number_of_elements = x_param.number_of_elements_;
      const std::size_t cuda_matrix_vector_length = ALU_param.vector_length_;
      const std::size_t number_of_groups = (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;

      // Shift the index for different groups
      d_ALU = d_ALU + group_id * ALU_param.number_of_elements_ / number_of_groups;
      d_x = d_x + group_id * x_param.number_of_elements_ / number_of_groups;
      double* d_y = d_x;  // Alias d_x for consistency with equation, but to reuse memory

      if (tid < number_of_grid_cells)
      {
        // Forward Substitution
        {
          for (auto i = 0; i < d_nLij_size; ++i)
          {
            const std::size_t j_lim = d_nLij[i];
            for (auto j = 0; j < j_lim; ++j)
            {
              const std::size_t d_Lij_yj_first = (*d_Lij_yj).first;
              const std::size_t d_Lij_yj_second_times_vector_length = (*d_Lij_yj).second * cuda_matrix_vector_length;
              auto d_ALU_ptr = d_ALU + d_Lij_yj_first;
              auto d_x_ptr = d_x + d_Lij_yj_second_times_vector_length;
              d_y[local_tid] -= d_ALU_ptr[local_tid] * d_x_ptr[local_tid];
              ++d_Lij_yj;
            }
            d_y += cuda_matrix_vector_length;
          }
        }
        // Backward Substitution
        {
          // d_y will be x_elem in the CPU implementation
          d_y = d_x + x_param.number_of_elements_ / number_of_groups - cuda_matrix_vector_length;
          for (auto i = 0; i < d_nUij_Uii_size; ++i)
          {
            const std::size_t j_lim = d_nUij_Uii[i].first;
            for (auto j = 0; j < j_lim; ++j)
            {
              auto d_ALU_ptr = d_ALU + (*d_Uij_xj).first;
              auto d_x_ptr = d_x + (*d_Uij_xj).second * cuda_matrix_vector_length;
              d_y[local_tid] -= d_ALU_ptr[local_tid] * d_x_ptr[local_tid];
              ++d_Uij_xj;
            }
            auto d_ALU_ptr = d_ALU + d_nUij_Uii[i].second;
            d_y[local_tid] /= d_ALU_ptr[local_tid];
            d_y -= cuda_matrix_vector_length;
          }
        }
      }
    }

    /// This is the function that will copy the constant data
    ///   members of class "CudaLinearSolver" to the device
    LinearSolverInPlaceParam CopyConstData(LinearSolverInPlaceParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t nLij_bytes = sizeof(size_t) * hoststruct.nLij_size_;
      size_t Lij_yj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.Lij_yj_size_;
      size_t nUij_Uii_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.nUij_Uii_size_;
      size_t Uij_xj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.Uij_xj_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      LinearSolverInPlaceParam devstruct;
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.nLij_), nLij_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.Lij_yj_), Lij_yj_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.nUij_Uii_), nUij_Uii_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.Uij_xj_), Uij_xj_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.nLij_,
              hoststruct.nLij_,
              nLij_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.Lij_yj_,
              hoststruct.Lij_yj_,
              Lij_yj_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.nUij_Uii_,
              hoststruct.nUij_Uii_,
              nUij_Uii_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.Uij_xj_,
              hoststruct.Uij_xj_,
              Uij_xj_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");

      devstruct.nLij_size_ = hoststruct.nLij_size_;
      devstruct.Lij_yj_size_ = hoststruct.Lij_yj_size_;
      devstruct.nUij_Uii_size_ = hoststruct.nUij_Uii_size_;
      devstruct.Uij_xj_size_ = hoststruct.Uij_xj_size_;
      devstruct.number_of_non_zeros_ = hoststruct.number_of_non_zeros_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaLinearSolver" on the device
    void FreeConstData(LinearSolverInPlaceParam& devstruct)
    {
      if (devstruct.nLij_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.nLij_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.Lij_yj_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.Lij_yj_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.nUij_Uii_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.nUij_Uii_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.Uij_xj_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.Uij_xj_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
    }

    void
    SolveKernelDriver(CudaMatrixParam& x_param, const CudaMatrixParam& ALU_param, const LinearSolverInPlaceParam& devstruct)
    {
      std::size_t number_of_blocks = (x_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SolveKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          x_param, ALU_param, devstruct);
    }
  }  // namespace cuda
}  // namespace micm