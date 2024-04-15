// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/cuda_matrix.cuh>
#include <micm/util/cuda_param.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <type_traits>
#include <cuda_runtime.h>

namespace micm
{
  template<class T, class OrderingPolicy>
  class CudaSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   private:
    CudaMatrixParam param_;

   public:
    CudaSparseMatrix()
        : SparseMatrix<T, OrderingPolicy>()
    {
      this->param_.d_data_ = nullptr;
    }

    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(std::is_same_v<T, double>)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
    }
    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
      this->param_.d_data_ = nullptr;
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(
        std::is_same_v<T, double>)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      return *this;
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      this->param_.d_data_ = nullptr;
      return *this;
    }

    CudaSparseMatrix(const CudaSparseMatrix& other) requires(std::is_same_v<T, double>)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), __FILE__, __LINE__, "cudaMemcpyDeviceToDevice");
    }

    CudaSparseMatrix(const CudaSparseMatrix& other)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_.d_data_ = nullptr;
    }

    CudaSparseMatrix(CudaSparseMatrix&& other) noexcept
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_.d_data_ = nullptr;
      std::swap(this->param_, other.param_);
    }

    CudaSparseMatrix& operator=(const CudaSparseMatrix& other)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(other);
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), __FILE__, __LINE__, "cudaMemcpyDeviceToDevice");
      return *this;
    }

    CudaSparseMatrix& operator=(CudaSparseMatrix&& other) noexcept
    {
      if (this != &other)
      {
        SparseMatrix<T, OrderingPolicy>::operator=(std::move(other));
        std::swap(this->param_, other.param_);
      }
      return *this;
    }

    ~CudaSparseMatrix() requires(std::is_same_v<T, double>)
    {
      std::cout << "Freeing device memory at address: " << this->param_.d_data_ << std::endl;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), __FILE__, __LINE__, "cudaFree");
      this->param_.d_data_ = nullptr;
    }

    ~CudaSparseMatrix()
    {
      this->param_.d_data_ = nullptr;
    }

    void CopyToDevice()
    {
       micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDevice(this->param_, this->data_), __FILE__, __LINE__, "cudaMemcpyHostToDevice");
    }

    void CopyToHost()
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToHost(this->param_, this->data_), __FILE__, __LINE__, "cudaMemcpyDeviceToHost");
    }

    CudaMatrixParam AsDeviceParam()
    {
      return this->param_;
    }
  };
}  // namespace micm
