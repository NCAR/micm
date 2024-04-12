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
      micm::cuda::MallocVector(this->param_, this->data_.size());
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
      micm::cuda::MallocVector(this->param_, this->data_.size());
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
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
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
    }

    CudaSparseMatrix& operator=(const CudaSparseMatrix& other)
    {
      return *this = CudaSparseMatrix(other);
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
      micm::cuda::FreeVector(this->param_);
    }

    ~CudaSparseMatrix()
    {
    }

    int CopyToDevice()
    {
      return micm::cuda::CopyToDevice(this->param_, this->data_);
    }
    int CopyToHost()
    {
      return micm::cuda::CopyToHost(this->param_, this->data_);
    }
    CudaMatrixParam AsDeviceParam()
    {
      return this->param_;
    }
  };
}  // namespace micm
