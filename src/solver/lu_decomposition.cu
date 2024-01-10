// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <chrono>
#include <iostream>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs LU decomposition on the device
    /// Note that passing the reference "LuDecomposeParam&" will pass the
    ///   compilation but the execution of this CUDA test hangs somehow
    __global__ void DecomposeKernel(const double* d_A, double* d_L, double* d_U, LuDecomposeParam devstruct, size_t ngrids)
    {
      /// Local device variables
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

      std::pair<size_t, size_t>* d_niLU = devstruct.niLU_;
      char* d_do_aik = devstruct.do_aik_;
      size_t* d_aik = devstruct.aik_;
      std::pair<size_t, size_t>* d_uik_nkj = devstruct.uik_nkj_;
      std::pair<size_t, size_t>* d_lij_ujk = devstruct.lij_ujk_;
      char* d_do_aki = devstruct.do_aki_;
      size_t* d_aki = devstruct.aki_;
      std::pair<size_t, size_t>* d_lki_nkj = devstruct.lki_nkj_;
      std::pair<size_t, size_t>* d_lkj_uji = devstruct.lkj_uji_;
      size_t* d_uii = devstruct.uii_;
      size_t niLU_size = devstruct.niLU_size_;

      size_t do_aik_offset = 0;
      size_t aik_offset = 0;
      size_t uik_nkj_offset = 0;
      size_t lij_ujk_offset = 0;
      size_t do_aki_offset = 0;
      size_t aki_offset = 0;
      size_t lkj_uji_offset = 0;
      size_t lki_nkj_offset = 0;
      size_t uii_offset = 0;

      if (tid < ngrids)
      {
        // loop through every element in niLU
        for (size_t i = 0; i < niLU_size; i++)
        {
          // upper triangular matrix
          auto inLU = d_niLU[i];
          for (size_t iU = 0; iU < inLU.second; ++iU)
          {
            if (d_do_aik[do_aik_offset++])
            {
              size_t U_idx = d_uik_nkj[uik_nkj_offset].first + tid;
              size_t A_idx = d_aik[aik_offset++] + tid;
              d_U[U_idx] = d_A[A_idx];
            }

            for (size_t ikj = 0; ikj < d_uik_nkj[uik_nkj_offset].second; ++ikj)
            {
              size_t U_idx_1 = d_uik_nkj[uik_nkj_offset].first + tid;
              size_t L_idx = d_lij_ujk[lij_ujk_offset].first + tid;
              size_t U_idx_2 = d_lij_ujk[lij_ujk_offset].second + tid;
              d_U[U_idx_1] -= d_L[L_idx] * d_U[U_idx_2];
              ++lij_ujk_offset;
            }
            ++uik_nkj_offset;
          }
          // lower triangular matrix

          d_L[d_lki_nkj[lki_nkj_offset++].first + tid] = 1.0;

          for (size_t iL = 0; iL < inLU.first; ++iL)
          {
            if (d_do_aki[do_aki_offset++])
            {
              size_t L_idx = d_lki_nkj[lki_nkj_offset].first + tid;
              size_t A_idx = d_aki[aki_offset++] + tid;
              d_L[L_idx] = d_A[A_idx];
            }
            for (size_t ikj = 0; ikj < d_lki_nkj[lki_nkj_offset].second; ++ikj)
            {
              size_t L_idx_1 = d_lki_nkj[lki_nkj_offset].first + tid;
              size_t L_idx_2 = d_lkj_uji[lkj_uji_offset].first + tid;
              size_t U_idx = d_lkj_uji[lkj_uji_offset].second + tid;
              d_L[L_idx_1] -= d_L[L_idx_2] * d_U[U_idx];
              ++lkj_uji_offset;
            }
            d_L[d_lki_nkj[lki_nkj_offset].first + tid] /= d_U[d_uii[uii_offset] + tid];
            ++lki_nkj_offset;
            ++uii_offset;
          }
        }
      }
    }  // end of CUDA kernel

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecomposition" to the device
    LuDecomposeParam CopyConstData(LuDecomposeParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      size_t niLU_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.niLU_size_;
      size_t do_aik_bytes = sizeof(char) * hoststruct.do_aik_size_;
      size_t aik_bytes = sizeof(size_t) * hoststruct.aik_size_;
      size_t uik_nkj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.uik_nkj_size_;
      size_t lij_ujk_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.lij_ujk_size_;
      size_t do_aki_bytes = sizeof(char) * hoststruct.do_aki_size_;
      size_t aki_bytes = sizeof(size_t) * hoststruct.aki_size_;
      size_t lki_nkj_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.lki_nkj_size_;
      size_t lkj_uji_bytes = sizeof(std::pair<size_t, size_t>) * hoststruct.lkj_uji_size_;
      size_t uii_bytes = sizeof(size_t) * hoststruct.uii_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      LuDecomposeParam devstruct;
      cudaMalloc(&(devstruct.niLU_), niLU_bytes);
      cudaMalloc(&(devstruct.do_aik_), do_aik_bytes);
      cudaMalloc(&(devstruct.aik_), aik_bytes);
      cudaMalloc(&(devstruct.uik_nkj_), uik_nkj_bytes);
      cudaMalloc(&(devstruct.lij_ujk_), lij_ujk_bytes);
      cudaMalloc(&(devstruct.do_aki_), do_aki_bytes);
      cudaMalloc(&(devstruct.aki_), aki_bytes);
      cudaMalloc(&(devstruct.lki_nkj_), lki_nkj_bytes);
      cudaMalloc(&(devstruct.lkj_uji_), lkj_uji_bytes);
      cudaMalloc(&(devstruct.uii_), uii_bytes);

      /// Copy the data from host to device
      cudaMemcpy(devstruct.niLU_, hoststruct.niLU_, niLU_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.do_aik_, hoststruct.do_aik_, do_aik_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.aik_, hoststruct.aik_, aik_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.uik_nkj_, hoststruct.uik_nkj_, uik_nkj_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.lij_ujk_, hoststruct.lij_ujk_, lij_ujk_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.do_aki_, hoststruct.do_aki_, do_aki_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.aki_, hoststruct.aki_, aki_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.lki_nkj_, hoststruct.lki_nkj_, lki_nkj_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.lkj_uji_, hoststruct.lkj_uji_, lkj_uji_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devstruct.uii_, hoststruct.uii_, uii_bytes, cudaMemcpyHostToDevice);
      devstruct.niLU_size_ = hoststruct.niLU_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeParam& devstruct)
    {
      cudaFree(devstruct.niLU_);
      cudaFree(devstruct.do_aik_);
      cudaFree(devstruct.aik_);
      cudaFree(devstruct.uik_nkj_);
      cudaFree(devstruct.lij_ujk_);
      cudaFree(devstruct.do_aki_);
      cudaFree(devstruct.aki_);
      cudaFree(devstruct.lki_nkj_);
      cudaFree(devstruct.lkj_uji_);
      cudaFree(devstruct.uii_);
    }

    std::chrono::nanoseconds DecomposeKernelDriver(CudaSparseMatrixParam& sparseMatrix, const LuDecomposeParam& devstruct)
    {
      /// Create device pointers
      double* d_A;
      double* d_L;
      double* d_U;

      /// Allocate device memory
      cudaMalloc(&d_A, sizeof(double) * sparseMatrix.A_size_);
      cudaMalloc(&d_L, sizeof(double) * sparseMatrix.L_size_);
      cudaMalloc(&d_U, sizeof(double) * sparseMatrix.U_size_);

      /// Copy data from host to device
      cudaMemcpy(d_A, sparseMatrix.A_, sizeof(double) * sparseMatrix.A_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_L, sparseMatrix.L_, sizeof(double) * sparseMatrix.L_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_U, sparseMatrix.U_, sizeof(double) * sparseMatrix.U_size_, cudaMemcpyHostToDevice);

      size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      /// Call CUDA kernel and measure the execution time
      auto startTime = std::chrono::high_resolution_clock::now();
      DecomposeKernel<<<num_block, BLOCK_SIZE>>>(d_A, d_L, d_U, devstruct, sparseMatrix.n_grids_);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      /// Copy the data from device to host
      cudaMemcpy(sparseMatrix.L_, d_L, sizeof(double) * sparseMatrix.L_size_, cudaMemcpyDeviceToHost);
      cudaMemcpy(sparseMatrix.U_, d_U, sizeof(double) * sparseMatrix.U_size_, cudaMemcpyDeviceToHost);

      /// Clean up
      cudaFree(d_A);
      cudaFree(d_L);
      cudaFree(d_U);

      return kernel_duration;
    }  // end of DecomposeKernelDriver
  }    // end of namespace cuda
}  // end of namespace micm
