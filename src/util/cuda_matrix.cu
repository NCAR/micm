// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/util/types.hpp>

#include <cuda_runtime.h>

#include <vector>

namespace micm::cuda
{
  template<typename T>
  cudaError_t MallocArray(T*& array, Index num_elements)
  {
    cudaError_t err =
        cudaMallocAsync(&array, sizeof(T) * num_elements, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    return err;
  }

  template<typename T>
  cudaError_t MallocVector(CudaMatrixParam& param, Index number_of_elements)
  {
    param.number_of_elements_ = number_of_elements;
    cudaError_t err = cudaMallocAsync(
        &(param.d_data_), sizeof(T) * number_of_elements, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    return err;
  }

  template<typename T>
  cudaError_t FreeArray(T*& array)
  {
    if (array == nullptr)
    {
      return cudaError_t::cudaSuccess;
    }
    cudaError_t err = cudaFreeAsync(array, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    array = nullptr;
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
    cudaError_t err = cudaFreeAsync(param.d_data_, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    param.d_data_ = nullptr;
    return err;
  }

  template<typename T>
  __global__ void MatrixMaxKernel(T* d_data, Index number_of_elements, T val)
  {
    Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if (tid < number_of_elements)
    {
      d_data[tid] = d_data[tid] > val ? d_data[tid] : val;
    }
  }

  template<typename T>
  cudaError_t MatrixMax(CudaMatrixParam& param, T val)
  {
    const Index number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
    MatrixMaxKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
        param.d_data_, param.number_of_elements_, val);
    cudaError_t err = cudaGetLastError();
    return err;
  }

  template<typename T>
  __global__ void MatrixMinKernel(T* d_data, Index number_of_elements, T val)
  {
    Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if (tid < number_of_elements)
    {
      d_data[tid] = d_data[tid] < val ? d_data[tid] : val;
    }
  }

  template<typename T>
  cudaError_t MatrixMin(CudaMatrixParam& param, T val)
  {
    const Index number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
    MatrixMinKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
        param.d_data_, param.number_of_elements_, val);
    cudaError_t err = cudaGetLastError();
    return err;
  }

  template<typename T>
  cudaError_t CopyToDevice(CudaMatrixParam& param, const std::vector<T>& h_data)
  {
    cudaError_t err = cudaMemcpyAsync(
        param.d_data_,
        h_data.data(),
        sizeof(T) * param.number_of_elements_,
        cudaMemcpyHostToDevice,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    return err;
  }

  template<typename T>
  cudaError_t CopyToHost(CudaMatrixParam& param, std::vector<T>& h_data)
  {
    cudaError_t err = cudaMemcpyAsync(
        h_data.data(),
        param.d_data_,
        sizeof(T) * param.number_of_elements_,
        cudaMemcpyDeviceToHost,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    cudaStreamSynchronize(micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    return err;
  }

  template<typename T>
  cudaError_t CopyToDeviceFromDevice(CudaMatrixParam& vectorMatrixDest, const CudaMatrixParam& vectorMatrixSrc)
  {
    cudaError_t err = cudaMemcpyAsync(
        vectorMatrixDest.d_data_,
        vectorMatrixSrc.d_data_,
        sizeof(T) * vectorMatrixSrc.number_of_elements_,
        cudaMemcpyDeviceToDevice,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
    return err;
  }

  template<typename T>
  __global__ void FillCudaMatrixKernel(T* d_data, Index number_of_elements, T val)
  {
    Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if (tid < number_of_elements)
    {
      d_data[tid] = val;
    }
  }

  template<typename T>
  cudaError_t FillCudaMatrix(CudaMatrixParam& param, T val)
  {
    const Index number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
    FillCudaMatrixKernel<<<
        number_of_blocks,
        BLOCK_SIZE,
        0,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(param.d_data_, param.number_of_elements_, val);
    cudaError_t err = cudaGetLastError();
    return err;
  }

  // source code needs the instantiation of the template
  template cudaError_t MallocArray<Real>(Real*& array, Index num_elements);
  template cudaError_t MallocArray<Index>(Index*& array, Index num_elements);
  template cudaError_t FreeArray<Real>(Real*& array);
  template cudaError_t FreeArray<Index>(Index*& array);
  template cudaError_t MallocVector<Real>(CudaMatrixParam& param, Index number_of_elements);
  template cudaError_t MallocVector<int>(CudaMatrixParam& param, Index number_of_elements);
  template cudaError_t MatrixMax<Real>(CudaMatrixParam& param, Real val);
  template cudaError_t MatrixMin<Real>(CudaMatrixParam& param, Real val);
  template cudaError_t CopyToDevice<Real>(CudaMatrixParam& param, const std::vector<Real>& h_data);
  template cudaError_t CopyToHost<Real>(CudaMatrixParam& param, std::vector<Real>& h_data);
  template cudaError_t CopyToDeviceFromDevice<Real>(
      CudaMatrixParam& vectorMatrixDest,
      const CudaMatrixParam& vectorMatrixSrc);
  template cudaError_t CopyToDeviceFromDevice<int>(
      CudaMatrixParam& vectorMatrixDest,
      const CudaMatrixParam& vectorMatrixSrc);
  template cudaError_t FillCudaMatrix<Real>(CudaMatrixParam& param, Real val);
}  // namespace micm::cuda
