// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <chrono>
#include <iostream>
#include <micm/util/cuda_param.hpp>
#include <vector>
#include <micm/solver/cuda_lu_decomposition.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs LU decomposition on the device
    __global__ void DecomposeKernel(const double* d_A, double* d_L, double* d_U, 
                                    LuDecomposeConstDevice* devptr, size_t ngrids)
    {
      /// Local device variables
      size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      std::pair<size_t, size_t>* uik_nkj = devptr->d_uik_nkj_;
      std::pair<size_t, size_t>* lij_ujk = devptr->d_lij_ujk_;
      std::pair<size_t, size_t>* lkj_uji = devptr->d_lkj_uji_;
      std::pair<size_t, size_t>* lki_nkj = device->d_lki_nkj_;
      size_t niLu_size = devptr->d_niLU_.size();
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
          auto inLU = device->d_niLU_[i];
          for (size_t iU = 0; iU < inLU.second; ++iU)
          {
            if (device->d_do_aik_[do_aik_offset++])
            {
              size_t U_idx = uik_nkj[uik_nkj_offset].first + tid;
              size_t A_idx = device->d_aik_[aik_offset++] + tid;
              d_U[U_idx] = d_A[A_idx];
            }

            for (size_t ikj = 0; ikj < uik_nkj[uik_nkj_offset].second; ++ikj)
            {
              size_t U_idx_1 = uik_nkj[uik_nkj_offset].first + tid;
              size_t L_idx = lij_ujk[lij_ujk_offset].first + tid;
              size_t U_idx_2 = lij_ujk[lij_ujk_offset].second + tid;
              d_U[U_idx_1] -= d_L[L_idx] * d_U[U_idx_2];
              ++lij_ujk_offset;
            }
            ++uik_nkj_offset;
          }
          // lower triangular matrix

          d_L[lki_nkj[lki_nkj_offset++].first + tid] = 1.0;

          for (size_t iL = 0; iL < inLU.first; ++iL)
          {
            if (device->d_do_aki_[do_aki_offset++])
            {
              size_t L_idx = lki_nkj[lki_nkj_offset].first + tid;
              size_t A_idx = device->d_aki_[aki_offset++] + tid;
              d_L[L_idx] = d_A[A_idx];
            }
            for (size_t ikj = 0; ikj < lki_nkj[lki_nkj_offset].second; ++ikj)
            {
              size_t L_idx_1 = lki_nkj[lki_nkj_offset].first + tid;
              size_t L_idx_2 = lkj_uji[lkj_uji_offset].first + tid;
              size_t U_idx = lkj_uji[lkj_uji_offset].second + tid;
              d_L[L_idx_1] -= d_L[L_idx_2] * d_U[U_idx];
              ++lkj_uji_offset;
            }
            d_L[lki_nkj[lki_nkj_offset].first + tid] /= d_U[device->uii_[uii_offset] + tid];
            ++lki_nkj_offset;
            ++uii_offset;
          }
        }
      }
    }  // end of CUDA kernel

    /// This is the function that copies the constant data members
    /// of objects with the "CudaLuDecomposition" type to the device
    void CopyConstData(CudaLuDecomposition* self, LuDecomposeConstDevice* devptr)
    {
      /// allocate device memory for the device struct
      cudaMalloc(&devptr,               sizeof(LuDecomposeConstDevice));
      cudaMalloc(&(devptr->d_niLU_),    sizeof(std::pair<size_t, size_t>) * self.niLU_.size());
      cudaMalloc(&(devptr->d_do_aik_),  sizeof(char) * self.do_aik_.size());
      cudaMalloc(&(devptr->d_aik_),     sizeof(size_t) * self.aik_.size());
      cudaMalloc(&(devptr->d_uik_nkj_), sizeof(std::pair<size_t, size_t>) * self.uik_nkj_.size());
      cudaMalloc(&(devptr->d_lij_ujk_), sizeof(std::pair<size_t, size_t>) * self.lij_ujk_.size());
      cudaMalloc(&(devptr->d_do_aki_),  sizeof(char) * self.do_aki_.size());
      cudaMalloc(&(devptr->d_aki_),     sizeof(size_t) * self.aki_.size());
      cudaMalloc(&(devptr->d_lki_nkj_), sizeof(std::pair<size_t, size_t>) * self.lki_nkj_.size());
      cudaMalloc(&(devptr->d_lkj_uji_), sizeof(std::pair<size_t, size_t>) * self.lkj_uji_.size());
      cudaMalloc(&(devptr->d_uii_),     sizeof(size_t) * self.uii_.size());

      /// copy the data from host to device
      cudaMemcpy(&(devptr->d_niLU_),    self.niLU_.data(), sizeof(std::pair<size_t, size_t>) * self.niLU_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_do_aik_),  self.do_aik_.data(), sizeof(char) * self.do_aik_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_aik_),     self.aik_.data(), sizeof(size_t) * self.aik_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_uik_nkj_), self.uik_nkj_.data(), sizeof(std::pair<size_t, size_t>) * self.uik_nkj_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_lij_ujk_), self.lij_ujk_.data(), sizeof(std::pair<size_t, size_t>) * self.lij_ujk_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_do_aki_),  self.do_aki_.data(), sizeof(char) * self.do_aki_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_aki_),     self.aki_.data(), sizeof(size_t) * self.aki_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_lki_nkj_), self.lki_nkj_.data(), sizeof(std::pair<size_t, size_t>) * self.lki_nkj_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_lkj_uji_), self.lkj_uji_.data(), sizeof(std::pair<size_t, size_t>) * self.lkj_uji_.size(), cudaMemcpyHostToDevice);
      cudaMemcpy(&(devptr->d_uii_),     self.uii_.data(), sizeof(size_t) * self.uii_.size(), cudaMemcpyHostToDevice);
    }

    /// This is the function that deletes the constant data members
    /// of objects with the "CudaLuDecomposition" type on the device
    void FreeConstData(LuDecomposeConstDevice* devptr)
    {
      cudaFree(devptr);
    }

    std::chrono::nanoseconds DecomposeKernelDriver(CudaSparseMatrixParam& sparseMatrix, LuDecomposeConstDevice* devptr)
    {
      /// create device pointers and allocate device memory
      double* d_A;
      double* d_L;
      double* d_U;

      cudaMalloc(&d_A, sizeof(double) * sparseMatrix.A_size_);
      cudaMalloc(&d_L, sizeof(double) * sparseMatrix.L_size_);
      cudaMalloc(&d_U, sizeof(double) * sparseMatrix.U_size_);

      /// copy data from host to device
      cudaMemcpy(d_A, sparseMatrix.A_, sizeof(double) * sparseMatrix.A_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_L, sparseMatrix.L_, sizeof(double) * sparseMatrix.L_size_, cudaMemcpyHostToDevice);
      cudaMemcpy(d_U, sparseMatrix.U_, sizeof(double) * sparseMatrix.U_size_, cudaMemcpyHostToDevice);

      size_t num_block = (sparseMatrix.n_grids_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      /// call CUDA kernel and measure the execution time
      auto startTime = std::chrono::high_resolution_clock::now();
      DecomposeKernel<<<num_block, BLOCK_SIZE>>>(d_A, d_L, d_U, devptr, sparseMatrix.n_grids_);
      cudaDeviceSynchronize();
      auto endTime = std::chrono::high_resolution_clock::now();
      auto kernel_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

      /// copy the data from device to host
      cudaMemcpy(sparseMatrix.L_, d_L, sizeof(double) * sparseMatrix.L_size_, cudaMemcpyDeviceToHost);
      cudaMemcpy(sparseMatrix.U_, d_U, sizeof(double) * sparseMatrix.U_size_, cudaMemcpyDeviceToHost);
      
      /// clean up
      cudaFree(d_A);
      cudaFree(d_L);
      cudaFree(d_U);

      return kernel_duration;
    }  // end of DecomposeKernelDriver
  }    // end of namespace cuda
}      // end of namespace micm
