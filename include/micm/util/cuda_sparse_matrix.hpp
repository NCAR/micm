// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/cuda_matrix.cuh>
#include <micm/util/cuda_param.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <type_traits>

namespace micm
{
  template<class T, class OrderingPolicy>
  class CudaSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   private:
    CudaMatrixParam param_;

   public:
    CudaSparseMatrix() = default;

    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(std::is_same_v<T, double>)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      cudaError_t err = micm::cuda::MallocVector(this->param_, this->data_.size());
      if (err != cudaSuccess)
      {
          throw std::runtime_error("cudaMalloc failed: " + std::string(cudaGetErrorString(err)));
      }
    }
    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(
        std::is_same_v<T, double>)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      cudaError_t err = micm::cuda::MallocVector(this->param_, this->data_.size());
      if (err != cudaSuccess)
      {
          throw std::runtime_error("cudaMalloc failed: " + std::string(cudaGetErrorString(err)));
      }
      return *this;
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      return *this;
    }

    CudaSparseMatrix(const CudaSparseMatrix& other) requires(std::is_same_v<T, double>)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::MallocVector(this->param_, this->data_.size());
      micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_);
    }

    CudaSparseMatrix(const CudaSparseMatrix& other)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
    }

    CudaSparseMatrix(CudaSparseMatrix&& other) noexcept
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_ = std::move(other.param_);
      other.param_.d_data_ = nullptr;
    }

    CudaSparseMatrix& operator=(const CudaSparseMatrix& other)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(other);
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::MallocVector(this->param_, this->data_.size());
      micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_); 
      return *this;
    }

    CudaSparseMatrix& operator=(CudaSparseMatrix&& other) noexcept
    {
      if (this != &other)
      {
        SparseMatrix<T, OrderingPolicy>::operator=(std::move(other));
        this->param_ = std::move(other.param_);
        other.param_.d_data_ = nullptr;
      }
      return *this;
    }

    ~CudaSparseMatrix() requires(std::is_same_v<T, double>)
    {
      cudaError_t err = micm::cuda::FreeVector(this->param_);
      if (err != cudaSuccess)
      {
          throw std::runtime_error("cudaFree failed: " + std::string(cudaGetErrorString(err)));
      }
    }

    ~CudaSparseMatrix()
    {
    }

    void CopyToDevice()
    {
      cudaError_t err = micm::cuda::CopyToDevice(this->param_, this->data_);
      if (err != cudaSuccess)
      {
          throw std::runtime_error("cudaMemcpyHostToDevice failed: " + std::string(cudaGetErrorString(err)));
      }
    }
    void CopyToHost()
    {
      cudaError_t err = micm::cuda::CopyToHost(this->param_, this->data_);
      if (err != cudaSuccess)
      {
          throw std::runtime_error("cudaMemcpyDeviceToHost failed: " + std::string(cudaGetErrorString(err)));
      }
    }
    CudaMatrixParam AsDeviceParam()
    {
      return this->param_;
    }
  };
}  // namespace micm
