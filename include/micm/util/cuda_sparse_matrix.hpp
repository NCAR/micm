/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/util/cuda_dense_matrix.hpp>  // include this for CudaMatrix concept
#include <micm/util/cuda_matrix.cuh>
#include <micm/util/cuda_param.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cuda_runtime.h>

#include <type_traits>

namespace micm
{
  template<class T, class OrderingPolicy>
  class CudaSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = CudaSparseMatrix<int, OrderingPolicy>;
    using value_type = T;

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
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), "cudaMalloc");
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
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), "cudaMalloc");
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
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
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
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
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
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), "cudaFree");
      this->param_.d_data_ = nullptr;
    }

    ~CudaSparseMatrix()
    {
      this->param_.d_data_ = nullptr;
    }

    void CopyToDevice()
    {
      static_assert(std::is_same_v<T, double>);
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice(this->param_, this->data_), "cudaMemcpyHostToDevice");
    }

    void CopyToHost()
    {
      static_assert(std::is_same_v<T, double>);
      CHECK_CUDA_ERROR(micm::cuda::CopyToHost(this->param_, this->data_), "cudaMemcpyDeviceToHost");
    }

    CudaMatrixParam AsDeviceParam() const
    {
      return this->param_;
    }
  };
}  // namespace micm
