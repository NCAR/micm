// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs LU decomposition on the device
    __global__ void DecomposeKernel(
        const CudaMatrixParam A_param,
        CudaMatrixParam L_param,
        CudaMatrixParam U_param,
        const LuDecomposeParam devstruct)
    {
      // Calculate global thread ID
      size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::pair<size_t, size_t>* const d_niLU = devstruct.niLU_;
      const char* const d_do_aik = devstruct.do_aik_;
      const size_t* const d_aik = devstruct.aik_;
      const std::pair<size_t, size_t>* const d_uik_nkj = devstruct.uik_nkj_;
      const std::pair<size_t, size_t>* const d_lij_ujk = devstruct.lij_ujk_;
      const char* const d_do_aki = devstruct.do_aki_;
      const size_t* const d_aki = devstruct.aki_;
      const std::pair<size_t, size_t>* const d_lki_nkj = devstruct.lki_nkj_;
      const std::pair<size_t, size_t>* const d_lkj_uji = devstruct.lkj_uji_;
      const size_t* const d_uii = devstruct.uii_;
      const size_t niLU_size = devstruct.niLU_size_;

      size_t do_aik_offset = 0;
      size_t aik_offset = 0;
      size_t uik_nkj_offset = 0;
      size_t lij_ujk_offset = 0;
      size_t do_aki_offset = 0;
      size_t aki_offset = 0;
      size_t lkj_uji_offset = 0;
      size_t lki_nkj_offset = 0;
      size_t uii_offset = 0;

      const double* const d_A = A_param.d_data_;
      double* const d_L = L_param.d_data_;
      double* const d_U = U_param.d_data_;
      const size_t number_of_grid_cells = A_param.number_of_grid_cells_;

      if (tid < number_of_grid_cells)
      {
        // loop through every element in niLU
        for (size_t i = 0; i < niLU_size; i++)
        {
          // upper triangular matrix
          auto inLU = d_niLU[i];
          for (size_t iU = 0; iU < inLU.second; ++iU)
          {
            size_t U_idx = d_uik_nkj[uik_nkj_offset].first + tid;
            if (d_do_aik[do_aik_offset++])
            {
              size_t A_idx = d_aik[aik_offset++] + tid;
              d_U[U_idx] = d_A[A_idx];
            }
            else
            {
              d_U[U_idx] = 0;
            }

            for (size_t ikj = 0; ikj < d_uik_nkj[uik_nkj_offset].second; ++ikj)
            {
              size_t L_idx = d_lij_ujk[lij_ujk_offset].first + tid;
              size_t U_idx_2 = d_lij_ujk[lij_ujk_offset].second + tid;
              d_U[U_idx] -= d_L[L_idx] * d_U[U_idx_2];
              ++lij_ujk_offset;
            }
            ++uik_nkj_offset;
          }
          // lower triangular matrix

          d_L[d_lki_nkj[lki_nkj_offset++].first + tid] = 1.0;

          for (size_t iL = 0; iL < inLU.first; ++iL)
          {
            size_t L_idx = d_lki_nkj[lki_nkj_offset].first + tid;
            if (d_do_aki[do_aki_offset++])
            {
              size_t A_idx = d_aki[aki_offset++] + tid;
              d_L[L_idx] = d_A[A_idx];
            }
            else
            {
              d_L[L_idx] = 0;
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
    ///   members of class "CudaDoolittleLuDecomposition" to the device
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
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.niLU_), niLU_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.do_aik_), do_aik_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.aik_), aik_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.uik_nkj_), uik_nkj_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.lij_ujk_), lij_ujk_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.do_aki_), do_aki_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.aki_), aki_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.lki_nkj_), lki_nkj_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.lkj_uji_), lkj_uji_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.uii_), uii_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.niLU_,
              hoststruct.niLU_,
              niLU_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.do_aik_,
              hoststruct.do_aik_,
              do_aik_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aik_,
              hoststruct.aik_,
              aik_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.uik_nkj_,
              hoststruct.uik_nkj_,
              uik_nkj_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.lij_ujk_,
              hoststruct.lij_ujk_,
              lij_ujk_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.do_aki_,
              hoststruct.do_aki_,
              do_aki_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aki_,
              hoststruct.aki_,
              aki_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.lki_nkj_,
              hoststruct.lki_nkj_,
              lki_nkj_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.lkj_uji_,
              hoststruct.lkj_uji_,
              lkj_uji_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.uii_,
              hoststruct.uii_,
              uii_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      devstruct.niLU_size_ = hoststruct.niLU_size_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaDoolittleLuDecomposition" on the device
    void FreeConstData(LuDecomposeParam& devstruct)
    {
      if (devstruct.niLU_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.niLU_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.do_aik_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.do_aik_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.aik_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.aik_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.uik_nkj_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.uik_nkj_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.lij_ujk_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.lij_ujk_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.do_aki_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.do_aki_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.aki_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.aki_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.lki_nkj_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.lki_nkj_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.lkj_uji_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.lkj_uji_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.uii_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.uii_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
    }

    void DecomposeKernelDriver(
        const CudaMatrixParam& A_param,
        CudaMatrixParam& L_param,
        CudaMatrixParam& U_param,
        const LuDecomposeParam& devstruct)
    {
      // Launch the CUDA kernel for LU decomposition
      size_t number_of_blocks = (A_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      DecomposeKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          A_param, L_param, U_param, devstruct);
    }  // end of DecomposeKernelDriver
  }    // end of namespace cuda
}  // end of namespace micm
