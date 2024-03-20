// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs the "solve" function on the device
    __global__ void SolveKernel(
        double* d_lower_matrix,
        double* d_upper_matrix,
        double* d_b,
        double* d_x,
        LinearSolverParam devstruct,
        size_t n_grids,
        size_t b_column_counts,
        size_t x_column_counts)
    {
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      double* d_y = d_x;  // Alias d_x for consistency with equation, but to reuse memory
      std::pair<size_t, size_t>* d_nLij_Lii = devstruct.nLij_Lii_;
      std::pair<size_t, size_t>* d_Lij_yj = devstruct.Lij_yj_;
      std::pair<size_t, size_t>* d_nUij_Uii = devstruct.nUij_Uii_;
      std::pair<size_t, size_t>* d_Uij_xj = devstruct.Uij_xj_;
      size_t nLij_Lii_size = devstruct.nLij_Lii_size_;
      size_t nUij_Uii_size = devstruct.nUij_Uii_size_;

      if (tid < n_grids)
      {
        size_t b_column_index = 0;
        size_t x_column_index = 0;
        size_t y_column_index = 0;
        size_t b_column_backward_index = b_column_counts - 1;
        size_t x_column_backward_index = x_column_counts - 1;
        size_t Lij_yj_index = 0;
        size_t Uij_xj_index = 0;

        for (size_t j = 0; j < nLij_Lii_size; ++j)
        {
          auto& nLij_Lii_element = d_nLij_Lii[j];
          d_y[y_column_index * n_grids + tid] = d_b[b_column_index++ * n_grids + tid];
          for (size_t i = 0; i < nLij_Lii_element.first; ++i)
          {
            size_t lower_matrix_index = d_Lij_yj[Lij_yj_index].first + tid;
            size_t y_index = d_Lij_yj[Lij_yj_index].second * n_grids + tid;
            d_y[y_column_index * n_grids + tid] -= d_lower_matrix[lower_matrix_index] * d_y[y_index];
            ++Lij_yj_index;
          }
          d_y[y_column_index++ * n_grids + tid] /= d_lower_matrix[nLij_Lii_element.second + tid];
        }

        for (size_t k = 0; k < nUij_Uii_size; ++k)
        {
          auto& nUij_Uii_element = d_nUij_Uii[k];

          for (size_t i = 0; i < nUij_Uii_element.first; ++i)
          {
            size_t upper_matrix_index = d_Uij_xj[Uij_xj_index].first + tid;
            size_t x_index = d_Uij_xj[Uij_xj_index].second * n_grids + tid;
            d_x[x_column_backward_index * n_grids + tid] -= d_upper_matrix[upper_matrix_index] * d_x[x_index];
            ++Uij_xj_index;
          }
          d_x[x_column_backward_index * n_grids + tid] /= d_upper_matrix[nUij_Uii_element.second + tid];

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
      cudaMalloc(&(devstruct.nLij_Lii_), nLij_Lii_bytes);
      cudaMalloc(&(devstruct.Lij_yj_), Lij_yj_bytes);
      cudaMalloc(&(devstruct.nUij_Uii_), nUij_Uii_bytes);
      cudaMalloc(&(devstruct.Uij_xj_), Uij_xj_bytes);

      /// Copy the data from host to device
      cudaMemcpy(devstruct.nLij_Lii_, hoststruct.nLij_Lii_, nLij_Lii_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.Lij_yj_, hoststruct.Lij_yj_, Lij_yj_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.nUij_Uii_, hoststruct.nUij_Uii_, nUij_Uii_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.Uij_xj_, hoststruct.Uij_xj_, Uij_xj_bytes, cudaMemcpyHostToDevice);

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
      cudaFree(devstruct.nLij_Lii_);
      cudaFree(devstruct.Lij_yj_);
      cudaFree(devstruct.nUij_Uii_);
      cudaFree(devstruct.Uij_xj_);
    }

    std::chrono::nanoseconds
    SolveKernelDriver(CudaSparseMatrixParam& sparseMatrix, CudaMatrixParam& denseMatrix, const LinearSolverParam& devstruct)
    {
      /// Create device pointers
      double* d_lower_matrix;
      double* d_upper_matrix;
      double* d_b;
      double* d_x;

      /// Allocate device memory
      cudaMalloc(&d_lower_matrix, sizeof(double) * sparseMatrix.lower_matrix_size_);
      cudaMalloc(&d_upper_matrix, sizeof(double) * sparseMatrix.upper_matrix_size_);
      cudaMalloc(&d_b, sizeof(double) * denseMatrix.b_size_);
      cudaMalloc(&d_x, sizeof(double) * denseMatrix.x_size_);

      /// Copy data from host to device
      cudaMemcpy(
          d_lower_matrix,
          sparseMatrix.lower_matrix_,
          sizeof(double) * sparseMatrix.lower_matrix_size_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(
          d_upper_matrix,
          sparseMatrix.upper_matrix_,
          sizeof(double) * sparseMatrix.upper_matrix_size_,
          cudaMemcpyHostToDevice);
      cudaMemcpy(d_b, denseMatrix.b_, sizeof(double) * denseMatrix.b_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_x, denseMatrix.x_, sizeof(double) * denseMatrix.x_size_, cudaMemcpyHostToDevice);

      size_t num_block = (denseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      /// Call CUDA kernel and measure the execution time
      auto startTime = std::chrono::high_resolution_clock::now();
      SolveKernel<<<num_block, BLOCK_SIZE>>>(
          d_lower_matrix,
          d_upper_matrix,
          d_b,
          d_x,
          devstruct,
          denseMatrix.n_grids_,
          denseMatrix.b_column_counts_,
          denseMatrix.x_column_counts_);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      /// Copy the data from device to host
      cudaMemcpy(denseMatrix.x_, d_x, sizeof(double) * denseMatrix.x_size_, cudaMemcpyDeviceToHost);

      /// Clean up
      cudaFree(d_lower_matrix);
      cudaFree(d_upper_matrix);
      cudaFree(d_b);
      cudaFree(d_x);

      return kernel_duration;
    }
  }  // namespace cuda
}  // namespace micm