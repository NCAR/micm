// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cuda_runtime.h>
#include <iostream>
#include <micm/util/cuda_matrix.cuh>
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
      cudaError_t err = cudaMemcpy(param.d_data_, h_data.data(), sizeof(double) * param.number_of_elements_, cudaMemcpyHostToDevice);
      return err;
    }
    cudaError_t CopyToHost(CudaMatrixParam& param, std::vector<double>& h_data)
    {
      cudaDeviceSynchronize();
      cudaError_t err = cudaMemcpy(h_data.data(), param.d_data_, sizeof(double) * param.number_of_elements_, cudaMemcpyDeviceToHost);
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
  }  // namespace cuda
}  // namespace micm
