// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/util/internal_error.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <cuda_runtime.h>

#include <vector>

namespace micm
{
  namespace cuda
  {
    cudaError_t MallocVector(CudaMatrixParam& param, std::size_t number_of_elements)
    {
      param.number_of_elements_ = number_of_elements;
      cudaError_t err = cudaMalloc(&(param.d_data_), sizeof(double) * number_of_elements);
      return err;
    }

    cudaError_t FreeVector(CudaMatrixParam& param)
    {
      param.number_of_elements_ = 0;
      param.number_of_grid_cells_ = 0;
      if (param.d_data_ == nullptr)
      {
        return cudaError_t::cudaSuccess;
      }
      cudaError_t err = cudaFree(param.d_data_);
      param.d_data_ = nullptr;
      return err;
    }

    cudaError_t CopyToDevice(CudaMatrixParam& param, std::vector<double>& h_data)
    {
      cudaError_t err =
          cudaMemcpy(param.d_data_, h_data.data(), sizeof(double) * param.number_of_elements_, cudaMemcpyHostToDevice);
      return err;
    }

    cudaError_t CopyToHost(CudaMatrixParam& param, std::vector<double>& h_data)
    {
      cudaDeviceSynchronize();
      cudaError_t err =
          cudaMemcpy(h_data.data(), param.d_data_, sizeof(double) * param.number_of_elements_, cudaMemcpyDeviceToHost);
      return err;
    }

    cudaError_t CopyToDeviceFromDevice(CudaMatrixParam& vectorMatrixDest, const CudaMatrixParam& vectorMatrixSrc)
    {
      cudaError_t err = cudaMemcpy(
          vectorMatrixDest.d_data_,
          vectorMatrixSrc.d_data_,
          sizeof(double) * vectorMatrixSrc.number_of_elements_,
          cudaMemcpyDeviceToDevice);
      return err;
    }

    template<typename T>
    __global__ void FillCudaMatrixKernel(T* d_data, std::size_t number_of_elements, T val)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid < number_of_elements)
      {
        d_data[tid] = val;
      }
    }

    template<typename T>
    cudaError_t FillCudaMatrix(CudaMatrixParam& param, T val)
    {
      std::size_t number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      FillCudaMatrixKernel<<<number_of_blocks, BLOCK_SIZE>>>(param.d_data_, param.number_of_elements_, val);
      cudaError_t err = cudaGetLastError();
      return err;
    }

    template cudaError_t FillCudaMatrix<double>(CudaMatrixParam& param, double val);
  }  // namespace cuda
}  // namespace micm
