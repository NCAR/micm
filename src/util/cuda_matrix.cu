// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/util/internal_error.hpp>

#include <cuda_runtime.h>

#include <vector>

namespace micm
{
  namespace cuda
  {
    template<typename T>
    cudaError_t MallocVector(CudaMatrixParam& param, std::size_t number_of_elements)
    {
      param.number_of_elements_ = number_of_elements;
      cudaError_t err = cudaMallocAsync(
          &(param.d_data_), sizeof(T) * number_of_elements, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0));
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
    __global__ void MatrixMaxKernel(T* d_data, std::size_t number_of_elements, T val)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid < number_of_elements)
      {
        d_data[tid] = d_data[tid] > val ? d_data[tid] : val;
      }
    }

    template<typename T>
    cudaError_t MatrixMax(CudaMatrixParam& param, T val)
    {
      std::size_t number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      MatrixMaxKernel<<<number_of_blocks, BLOCK_SIZE, 0, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          param.d_data_, param.number_of_elements_, val);
      cudaError_t err = cudaGetLastError();
      return err;
    }

    template<typename T>
    __global__ void MatrixMinKernel(T* d_data, std::size_t number_of_elements, T val)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      if (tid < number_of_elements)
      {
        d_data[tid] = d_data[tid] < val ? d_data[tid] : val;
      }
    }

    template<typename T>
    cudaError_t MatrixMin(CudaMatrixParam& param, T val)
    {
      std::size_t number_of_blocks = (param.number_of_elements_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
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
      FillCudaMatrixKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(param.d_data_, param.number_of_elements_, val);
      cudaError_t err = cudaGetLastError();
      return err;
    }

    // source code needs the instantiation of the template
    template cudaError_t MallocVector<double>(CudaMatrixParam& param, std::size_t number_of_elements);
    template cudaError_t MallocVector<int>(CudaMatrixParam& param, std::size_t number_of_elements);
    template cudaError_t MatrixMax<double>(CudaMatrixParam& param, double val);
    template cudaError_t MatrixMin<double>(CudaMatrixParam& param, double val);
    template cudaError_t CopyToDevice<double>(CudaMatrixParam& param, const std::vector<double>& h_data);
    template cudaError_t CopyToHost<double>(CudaMatrixParam& param, std::vector<double>& h_data);
    template cudaError_t CopyToDeviceFromDevice<double>(
        CudaMatrixParam& vectorMatrixDest,
        const CudaMatrixParam& vectorMatrixSrc);
    template cudaError_t CopyToDeviceFromDevice<int>(
        CudaMatrixParam& vectorMatrixDest,
        const CudaMatrixParam& vectorMatrixSrc);
    template cudaError_t FillCudaMatrix<double>(CudaMatrixParam& param, double val);
  }  // namespace cuda
}  // namespace micm
