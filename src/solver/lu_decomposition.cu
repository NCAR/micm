// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <chrono>
#include <iostream>
#include <micm/util/cuda_param.hpp>
#include <vector>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs LU decomposition on the device
    __global__ void DecomposeKernel(const double* d_A, double* d_L, double* d_U,
                                    std::pair<size_t, size_t>* d_niLU,
                                    char* d_do_aik, size_t* d_aik,
                                    std::pair<size_t, size_t>* d_uik_nkj,
                                    std::pair<size_t, size_t>* d_lij_ujk,
                                    char* d_do_aki, size_t* d_aki,
                                    std::pair<size_t, size_t>* d_lki_nkj,
                                    std::pair<size_t, size_t>* d_lkj_uji,
                                    size_t* d_uii,
                                    size_t niLU_size, size_t ngrids)
    {
      /// Local device variables
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
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
    void CopyConstData(LuDecomposeConst* hostptr, LuDecomposeConst*& devptr)
    {
      /// Calculate the memory space of each constant data member
      size_t niLU_bytes    = sizeof(std::pair<size_t, size_t>) * hostptr->niLU_size_;
      size_t do_aik_bytes  = sizeof(char) * hostptr->do_aik_size_;
      size_t aik_bytes     = sizeof(size_t) * hostptr->aik_size_; 
      size_t uik_nkj_bytes = sizeof(std::pair<size_t, size_t>) * hostptr->uik_nkj_size_; 
      size_t lij_ujk_bytes = sizeof(std::pair<size_t, size_t>) * hostptr->lij_ujk_size_;
      size_t do_aki_bytes  = sizeof(char) * hostptr->do_aki_size_;
      size_t aki_bytes     = sizeof(size_t) * hostptr->aki_size_;
      size_t lki_nkj_bytes = sizeof(std::pair<size_t, size_t>) * hostptr->lki_nkj_size_;
      size_t lkj_uji_bytes = sizeof(std::pair<size_t, size_t>) * hostptr->lkj_uji_size_;
      size_t uii_bytes     = sizeof(size_t) * hostptr->uii_size_;

      /// Can not use "cudaMalloc((void**)&devptr, sizeof(LuDecomposeConstDevice))" 
      ///   because host variable "devptr" will contain addresss in the device memory (which is ok), 
      ///   but "devptr->d_niLU_" becomes illegal since we can not access/deference the address 
      ///   in the device memory from the host code directly; see more discussion from:
      ///   https://forums.developer.nvidia.com/t/cudamalloc-and-structs-and-pointers-problem/12266/2

      /// The solution is to keep devptr and its members as host variables,
      ///   but its members contain the addresses in the device memory.
      devptr = new LuDecomposeConst;
      cudaMalloc(&(devptr->niLU_),     niLU_bytes);
      cudaMalloc(&(devptr->do_aik_),   do_aik_bytes);
      cudaMalloc(&(devptr->aik_),      aik_bytes);
      cudaMalloc(&(devptr->uik_nkj_),  uik_nkj_bytes);      
      cudaMalloc(&(devptr->lij_ujk_),  lij_ujk_bytes);
      cudaMalloc(&(devptr->do_aki_),   do_aki_bytes);
      cudaMalloc(&(devptr->aki_),      aki_bytes);
      cudaMalloc(&(devptr->lki_nkj_),  lki_nkj_bytes);
      cudaMalloc(&(devptr->lkj_uji_),  lkj_uji_bytes);
      cudaMalloc(&(devptr->uii_),      uii_bytes);

      /// Copy the data from host to device
      cudaMemcpy(devptr->niLU_,    hostptr->niLU_,    niLU_bytes,    cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->do_aik_,  hostptr->do_aik_,  do_aik_bytes,  cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->aik_,     hostptr->aik_,     aik_bytes,     cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->uik_nkj_, hostptr->uik_nkj_, uik_nkj_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->lij_ujk_, hostptr->lij_ujk_, lij_ujk_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->do_aki_,  hostptr->do_aki_,  do_aki_bytes,  cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->aki_,     hostptr->aki_,     aki_bytes,     cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->lki_nkj_, hostptr->lki_nkj_, lki_nkj_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->lkj_uji_, hostptr->lkj_uji_, lkj_uji_bytes, cudaMemcpyHostToDevice);
      cudaMemcpy(devptr->uii_,     hostptr->uii_,     uii_bytes,     cudaMemcpyHostToDevice);
      devptr->niLU_size_ = hostptr->niLU_size_;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeConst*& devptr)
    {
      cudaFree(devptr->niLU_);
      cudaFree(devptr->do_aik_);
      cudaFree(devptr->aik_);
      cudaFree(devptr->uik_nkj_);      
      cudaFree(devptr->lij_ujk_);
      cudaFree(devptr->do_aki_);
      cudaFree(devptr->aki_);
      cudaFree(devptr->lki_nkj_);
      cudaFree(devptr->lkj_uji_);
      cudaFree(devptr->uii_);
    }

    std::chrono::nanoseconds DecomposeKernelDriver(CudaSparseMatrixParam& sparseMatrix, 
                                                   LuDecomposeConst* devptr)
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
      DecomposeKernel<<<num_block, BLOCK_SIZE>>>(d_A, d_L, d_U,
                                                 devptr->niLU_,
                                                 devptr->do_aik_, devptr->aik_,
                                                 devptr->uik_nkj_, devptr->lij_ujk_,
                                                 devptr->do_aki_, devptr->aki_,
                                                 devptr->lki_nkj_, devptr->lkj_uji_,
                                                 devptr->uii_, devptr->niLU_size_, 
                                                 sparseMatrix.n_grids_);
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
}      // end of namespace micm