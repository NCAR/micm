// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#include <stdexcept>
#include <string>

namespace micm
{
  namespace cuda
  {
    // Maximum matrix dimension for constant memory (12KB for aii+nji+nki at n=1024)
    static const std::size_t MAX_N_CONST = 1024;

    // Constant memory for per-iteration index arrays (warp-uniform broadcast)
    __constant__ uint32_t c_aii[MAX_N_CONST];
    __constant__ uint32_t c_nji[MAX_N_CONST];
    __constant__ uint32_t c_nki[MAX_N_CONST];

    /// LU decomposition kernel: 1 thread per matrix, uint32_t SoA indices.
    /// Per-iteration arrays (aii, nji, nki) are in constant memory.
    /// Inner-loop index pairs are packed as uint2 for single 8-byte loads.
    __global__ void DecomposeKernel(CudaMatrixParam ALU_param, const LuDecomposeParam devstruct)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid >= ALU_param.number_of_grid_cells_)
        return;

      const uint32_t* __restrict__ d_aji = devstruct.aji_;
      // Cast packed uint32_t arrays to uint2* for single 8-byte loads
      const uint2* __restrict__ d_aik_njk = reinterpret_cast<const uint2*>(devstruct.aik_njk_packed_);
      const uint2* __restrict__ d_ajk_aji = reinterpret_cast<const uint2*>(devstruct.ajk_aji_packed_);

      const std::size_t vl = ALU_param.vector_length_;
      const std::size_t local_tid = tid % vl;
      const std::size_t group_id = tid / vl;
      double* d_ALU = ALU_param.d_data_ + group_id * static_cast<std::size_t>(devstruct.number_of_non_zeros_) * vl;

      uint32_t aji_pos = 0;
      uint32_t aik_pos = 0;
      uint32_t ajk_pos = 0;
      const uint32_t n = devstruct.n_;

      for (uint32_t i = 0; i < n; ++i)
      {
        // aii, nji, nki from constant memory (single-cycle broadcast to entire warp)
        const double Aii_inverse = 1.0 / d_ALU[c_aii[i] + local_tid];

        const uint32_t nji = c_nji[i];
#pragma unroll
        for (uint32_t ij = 0; ij < nji; ++ij)
        {
          d_ALU[d_aji[aji_pos] + local_tid] *= Aii_inverse;
          ++aji_pos;
        }

        const uint32_t nki = c_nki[i];
        for (uint32_t ik = 0; ik < nki; ++ik)
        {
          // Single 8-byte load: aik index + njk count
          const uint2 aik_njk = d_aik_njk[aik_pos];
          const double Aik = d_ALU[aik_njk.x + local_tid];
          const uint32_t njk = aik_njk.y;
#pragma unroll
          for (uint32_t ijk = 0; ijk < njk; ++ijk)
          {
            // Single 8-byte load: ajk index + aji_update index
            const uint2 ajk_aji = d_ajk_aji[ajk_pos];
            d_ALU[ajk_aji.x + local_tid] -= d_ALU[ajk_aji.y + local_tid] * Aik;
            ++ajk_pos;
          }
          ++aik_pos;
        }
      }
    }

    /// Copy constant data to device
    LuDecomposeMozartInPlaceParam CopyConstData(LuDecomposeMozartInPlaceParam& hoststruct)
    {
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      std::size_t aji_bytes = sizeof(uint32_t) * hoststruct.aji_size_;
      std::size_t aik_njk_packed_bytes = sizeof(uint32_t) * 2 * hoststruct.aik_njk_size_;
      std::size_t ajk_aji_packed_bytes = sizeof(uint32_t) * 2 * hoststruct.ajk_aji_size_;

      LuDecomposeParam devstruct;

      /// Copy per-iteration arrays (aii, nji, nki) to constant memory for warp broadcast
      if (hoststruct.n_ > MAX_N_CONST)
      {
        throw std::runtime_error(
            "LU decomposition matrix dimension n=" + std::to_string(hoststruct.n_) +
            " exceeds constant memory limit of " + std::to_string(MAX_N_CONST));
      }
      std::size_t aii_bytes = sizeof(uint32_t) * hoststruct.n_;
      std::size_t nji_bytes = sizeof(uint32_t) * hoststruct.n_;
      std::size_t nki_bytes = sizeof(uint32_t) * hoststruct.n_;
      CHECK_CUDA_ERROR(
          cudaMemcpyToSymbol(c_aii, hoststruct.aii_, aii_bytes, 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");
      CHECK_CUDA_ERROR(
          cudaMemcpyToSymbol(c_nji, hoststruct.nji_, nji_bytes, 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");
      CHECK_CUDA_ERROR(
          cudaMemcpyToSymbol(c_nki, hoststruct.nki_, nki_bytes, 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");

      /// Allocate and copy sequential arrays to global device memory
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.aji_), aji_bytes, stream), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.aik_njk_packed_), aik_njk_packed_bytes, stream), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.ajk_aji_packed_), ajk_aji_packed_bytes, stream), "cudaMalloc");

      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(devstruct.aji_, hoststruct.aji_, aji_bytes, cudaMemcpyHostToDevice, stream), "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aik_njk_packed_, hoststruct.aik_njk_packed_, aik_njk_packed_bytes, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.ajk_aji_packed_, hoststruct.ajk_aji_packed_, ajk_aji_packed_bytes, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");

      /// Copy the scalar data members
      devstruct.n_ = hoststruct.n_;
      devstruct.aji_size_ = hoststruct.aji_size_;
      devstruct.aik_njk_size_ = hoststruct.aik_njk_size_;
      devstruct.ajk_aji_size_ = hoststruct.ajk_aji_size_;
      devstruct.number_of_non_zeros_ = hoststruct.number_of_non_zeros_;

      return devstruct;
    }

    /// Free device memory
    void FreeConstData(LuDecomposeMozartInPlaceParam& devstruct)
    {
      // aii, nji, nki are in __constant__ memory — no dynamic free needed
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      if (devstruct.aji_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.aji_, stream), "cudaFree");
      if (devstruct.aik_njk_packed_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.aik_njk_packed_, stream), "cudaFree");
      if (devstruct.ajk_aji_packed_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.ajk_aji_packed_, stream), "cudaFree");
    }

    void DecomposeKernelDriver(CudaMatrixParam& ALU_param, const LuDecomposeParam& devstruct)
    {
      std::size_t num_cells = ALU_param.number_of_grid_cells_;
      std::size_t num_blocks = (num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      DecomposeKernel<<<num_blocks, BLOCK_SIZE, 0, stream>>>(ALU_param, devstruct);
    }  // end of DecomposeKernelDriver
  }  // end of namespace cuda
}  // end of namespace micm
