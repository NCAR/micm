// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/util/cuda_param.hpp>
#include <micm/util/cuda_util.cuh>

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
      std::pair<size_t, size_t>* d_nLij_Lii = devstruct.nLij_Lii_;
      std::pair<size_t, size_t>* d_Lij_yj = devstruct.Lij_yj_;
      std::pair<size_t, size_t>* d_nUij_Uii = devstruct.nUij_Uii_;
      std::pair<size_t, size_t>* d_Uij_xj = devstruct.Uij_xj_;
      size_t nLij_Lii_size = devstruct.nLij_Lii_size_;
      size_t nUij_Uii_size = devstruct.nUij_Uii_size_;

      double* d_L = L_param.d_data_;
      double* d_U = U_param.d_data_;
      double* d_x = x_param.d_data_;
      double* d_y = d_x;  // Alias d_x for consistency with equation, but to reuse memory
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
      CHECK_CUDA_ERROR(cudaMalloc(&(devstruct.nLij_Lii_), nLij_Lii_bytes), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMalloc(&(devstruct.Lij_yj_), Lij_yj_bytes), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMalloc(&(devstruct.nUij_Uii_), nUij_Uii_bytes), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMalloc(&(devstruct.Uij_xj_), Uij_xj_bytes), "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpy(devstruct.nLij_Lii_, hoststruct.nLij_Lii_, nLij_Lii_bytes, cudaMemcpyHostToDevice), "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpy(devstruct.Lij_yj_, hoststruct.Lij_yj_, Lij_yj_bytes, cudaMemcpyHostToDevice), "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpy(devstruct.nUij_Uii_, hoststruct.nUij_Uii_, nUij_Uii_bytes, cudaMemcpyHostToDevice), "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpy(devstruct.Uij_xj_, hoststruct.Uij_xj_, Uij_xj_bytes, cudaMemcpyHostToDevice), "cudaMemcpy");

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
      CHECK_CUDA_ERROR(cudaFree(devstruct.nLij_Lii_), "cudaFree");
      CHECK_CUDA_ERROR(cudaFree(devstruct.Lij_yj_), "cudaFree");
      CHECK_CUDA_ERROR(cudaFree(devstruct.nUij_Uii_), "cudaFree");
      CHECK_CUDA_ERROR(cudaFree(devstruct.Uij_xj_), "cudaFree");
    }

    void SolveKernelDriver(
        CudaMatrixParam& x_param,
        const CudaMatrixParam& L_param,
        const CudaMatrixParam& U_param,
        const LinearSolverParam& devstruct)
    {
      size_t number_of_blocks = (b_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SolveKernel<<<number_of_blocks, BLOCK_SIZE>>>(b_param, x_param, L_param, U_param, devstruct);
    }
  }  // namespace cuda
}  // namespace micm