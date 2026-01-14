// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that performs LU decomposition on the device
    __global__ void DecomposeKernel(CudaMatrixParam ALU_param, const LuDecomposeParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::tuple<std::size_t, std::size_t, std::size_t>* const __restrict__ d_aii_nji_nki = devstruct.aii_nji_nki_;
      const std::size_t* __restrict__ d_aji = devstruct.aji_;
      const std::pair<std::size_t, std::size_t>* __restrict__ d_aik_njk = devstruct.aik_njk_;
      const std::pair<std::size_t, std::size_t>* __restrict__ d_ajk_aji = devstruct.ajk_aji_;
      const std::size_t d_aii_nji_nki_size = devstruct.aii_nji_nki_size_;

      double* __restrict__ d_ALU = ALU_param.d_data_;
      const std::size_t number_of_grid_cells = ALU_param.number_of_grid_cells_;
      const std::size_t cuda_matrix_vector_length = ALU_param.vector_length_;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = std::floor(static_cast<double>(tid) / cuda_matrix_vector_length);

      // Shift the index for different groups
      d_ALU += group_id * devstruct.number_of_non_zeros_ * cuda_matrix_vector_length;

      if (tid < number_of_grid_cells)
      {
        for (int i = 0; i < d_aii_nji_nki_size; ++i)
        {
          auto& d_aii_nji_nki_elem = d_aii_nji_nki[i];
          auto d_Aii = d_ALU + std::get<0>(d_aii_nji_nki_elem);
          auto d_Aii_inverse = 1.0 / d_Aii[local_tid];
          for (int ij = 0; ij < std::get<1>(d_aii_nji_nki_elem); ++ij)
          {
            auto d_ALU_ji = d_ALU + *d_aji + local_tid;
            *d_ALU_ji *= d_Aii_inverse;
            ++d_aji;
          }
          for (int ik = 0; ik < std::get<2>(d_aii_nji_nki_elem); ++ik)
          {
            const std::size_t d_aik_njk_first = std::get<0>(*d_aik_njk);
            const std::size_t d_aik_njk_second = std::get<1>(*d_aik_njk);
            for (std::size_t ijk = 0; ijk < d_aik_njk_second; ++ijk)
            {
              auto d_ALU_first = d_ALU + d_ajk_aji->first + local_tid;
              auto d_ALU_second = d_ALU + d_ajk_aji->second + local_tid;
              auto d_ALU_aik = d_ALU + d_aik_njk_first + local_tid;
              *d_ALU_first -= *d_ALU_second * *d_ALU_aik;
              ++d_ajk_aji;
            }
            ++d_aik_njk;
          }
        }
      }
    }  // end of CUDA kernel

    /// This is the function that will copy the constant data
    ///   members of class "CudaDoolittleLuDecomposition" to the device
    LuDecomposeMozartInPlaceParam CopyConstData(LuDecomposeMozartInPlaceParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      std::size_t aii_nji_nki_bytes =
          sizeof(std::tuple<std::size_t, std::size_t, std::size_t>) * hoststruct.aii_nji_nki_size_;
      std::size_t aji_bytes = sizeof(std::size_t) * hoststruct.aji_size_;
      std::size_t aik_njk_bytes = sizeof(std::pair<std::size_t, std::size_t>) * hoststruct.aik_njk_size_;
      std::size_t ajk_aji_bytes = sizeof(std::pair<std::size_t, std::size_t>) * hoststruct.ajk_aji_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      LuDecomposeParam devstruct;
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.aii_nji_nki_), aii_nji_nki_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.aji_), aji_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.aik_njk_), aik_njk_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(devstruct.ajk_aji_), ajk_aji_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aii_nji_nki_,
              hoststruct.aii_nji_nki_,
              aii_nji_nki_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aji_,
              hoststruct.aji_,
              aji_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.aik_njk_,
              hoststruct.aik_njk_,
              aik_njk_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.ajk_aji_,
              hoststruct.ajk_aji_,
              ajk_aji_bytes,
              cudaMemcpyHostToDevice,
              micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMemcpy");

      /// Copy the other data members from host to device
      devstruct.aii_nji_nki_size_ = hoststruct.aii_nji_nki_size_;
      devstruct.aji_size_ = hoststruct.aji_size_;
      devstruct.aik_njk_size_ = hoststruct.aik_njk_size_;
      devstruct.ajk_aji_size_ = hoststruct.ajk_aji_size_;
      devstruct.number_of_non_zeros_ = hoststruct.number_of_non_zeros_;

      return devstruct;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaDoolittleLuDecomposition" on the device
    void FreeConstData(LuDecomposeMozartInPlaceParam& devstruct)
    {
      if (devstruct.aii_nji_nki_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.aii_nji_nki_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
            "cudaFree");
      if (devstruct.aji_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.aji_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.aik_njk_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.aik_njk_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
      if (devstruct.ajk_aji_ != nullptr)
        CHECK_CUDA_ERROR(
            cudaFreeAsync(devstruct.ajk_aji_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)), "cudaFree");
    }

    void DecomposeKernelDriver(CudaMatrixParam& ALU_param, const LuDecomposeParam& devstruct)
    {
      // Launch the CUDA kernel for LU decomposition
      const std::size_t number_of_blocks = (ALU_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      DecomposeKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          ALU_param, devstruct);
    }  // end of DecomposeKernelDriver
  }  // end of namespace cuda
}  // end of namespace micm
