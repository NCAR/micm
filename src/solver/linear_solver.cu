// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <chrono>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs the "solve" function on the device
    __global__ void SolveKernel(
        CudaMatrixParam x_param,
        const CudaMatrixParam L_param,
        const CudaMatrixParam U_param,
        const LinearSolverParam devstruct)
    {
      // Calculate global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::pair<size_t, size_t>* const d_nLij_Lii = devstruct.nLij_Lii_;
      const std::pair<size_t, size_t>* const d_Lij_yj = devstruct.Lij_yj_;
      const std::pair<size_t, size_t>* const d_nUij_Uii = devstruct.nUij_Uii_;
      const std::pair<size_t, size_t>* const d_Uij_xj = devstruct.Uij_xj_;
      const size_t nLij_Lii_size = devstruct.nLij_Lii_size_;
      const size_t nUij_Uii_size = devstruct.nUij_Uii_size_;

      const double* const d_L = L_param.d_data_;
      const double* const d_U = U_param.d_data_;
      double* const d_x = x_param.d_data_;
      double* const d_y = d_x;  // Alias d_x for consistency with equation, but to reuse memory
      const size_t number_of_grid_cells = x_param.number_of_grid_cells_;
      const size_t number_of_species = x_param.number_of_elements_ / number_of_grid_cells;

      if (tid < number_of_grid_cells)
      {
        size_t x_column_index = 0;
        size_t y_column_index = 0;
        size_t x_column_backward_index = number_of_species - 1;
        size_t Lij_yj_index = 0;
        size_t Uij_xj_index = 0;

        for (size_t j = 0; j < nLij_Lii_size; ++j)
        {
          auto& nLij_Lii_element = d_nLij_Lii[j];
          d_y[(y_column_index * number_of_grid_cells) + tid] = d_x[(x_column_index++ * number_of_grid_cells) + tid];
          for (size_t i = 0; i < nLij_Lii_element.first; ++i)
          {
            size_t lower_matrix_index = d_Lij_yj[Lij_yj_index].first + tid;
            size_t y_index = d_Lij_yj[Lij_yj_index].second * number_of_grid_cells + tid;
            d_y[(y_column_index * number_of_grid_cells) + tid] -= d_L[lower_matrix_index] * d_y[y_index];
            ++Lij_yj_index;
          }
          d_y[(y_column_index++ * number_of_grid_cells) + tid] /= d_L[nLij_Lii_element.second + tid];
        }

        for (size_t k = 0; k < nUij_Uii_size; ++k)
        {
          auto& nUij_Uii_element = d_nUij_Uii[k];

          for (size_t i = 0; i < nUij_Uii_element.first; ++i)
          {
            size_t upper_matrix_index = d_Uij_xj[Uij_xj_index].first + tid;
            size_t x_index = d_Uij_xj[Uij_xj_index].second * number_of_grid_cells + tid;
            d_x[(x_column_backward_index * number_of_grid_cells) + tid] -= d_U[upper_matrix_index] * d_x[x_index];
            ++Uij_xj_index;
          }
          d_x[(x_column_backward_index * number_of_grid_cells) + tid] /= d_U[nUij_Uii_element.second + tid];

          if (x_column_backward_index != 0)
          {
            --x_column_backward_index;
          }
        }
      }
    }

    /// This is the function that will copy the constant data
    ///   members of class "CudaLinearSolver" to the device
    LinearSolverParam CopyConstData(LinearSolverParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t nLij_Lii_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.nLij_Lii_size_;
      size_t Lij_yj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.Lij_yj_size_;
      size_t nUij_Uii_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.nUij_Uii_size_;
      size_t Uij_xj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.Uij_xj_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      LinearSolverParam devstruct;
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.nLij_Lii_), nLij_Lii_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
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
              devstruct.nLij_Lii_,
              hoststruct.nLij_Lii_,
              nLij_Lii_bytes,
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

      devstruct.nLij_Lii_size_ = hoststruct.nLij_Lii_size_;
      devstruct.Lij_yj_size_ = hoststruct.Lij_yj_size_;
      devstruct.nUij_Uii_size_ = hoststruct.nUij_Uii_size_;
      devstruct.Uij_xj_size_ = hoststruct.Uij_xj_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaLinearSolver" on the device
    void FreeConstData(LinearSolverParam& devstruct)
    {
      if (devstruct.nLij_Lii_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.nLij_Lii_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
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

    void SolveKernelDriver(
        CudaMatrixParam& x_param,
        const CudaMatrixParam& L_param,
        const CudaMatrixParam& U_param,
        const LinearSolverParam& devstruct)
    {
      size_t number_of_blocks = (x_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SolveKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          x_param, L_param, U_param, devstruct);
    }
  }  // namespace cuda
}  // namespace micm