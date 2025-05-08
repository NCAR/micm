// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>  // include this for CudaMatrix concept
#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>
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

    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector<T>(this->param_, this->data_.size()), "cudaMalloc");
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      this->param_.number_of_grid_cells_ = this->number_of_blocks_;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector<T>(this->param_, this->data_.size()), "cudaMalloc");
      return *this;
    }

    CudaSparseMatrix(const CudaSparseMatrix& other)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector<T>(this->param_, this->data_.size()), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice<T>(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
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
      CHECK_CUDA_ERROR(micm::cuda::MallocVector<T>(this->param_, this->data_.size()), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice<T>(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
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

    ~CudaSparseMatrix()
    {
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), "cudaFree");
      this->param_.d_data_ = nullptr;
    }

    void CopyToDevice()
    {
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice<T>(this->param_, this->data_), "cudaMemcpyHostToDevice");
    }

    void CopyToHost()
    {
      CHECK_CUDA_ERROR(micm::cuda::CopyToHost<T>(this->param_, this->data_), "cudaMemcpyDeviceToHost");
    }

    /// @brief Set every matrix element to a given value on the GPU
    /// @param val Value to set each element to
    void Fill(T val)
    {
      if constexpr (std::is_same_v<T, int>)
      {
        // the cudaMemset function only works for integer types and is an asynchronous function:
        // https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1gf7338650f7683c51ee26aadc6973c63a
        CHECK_CUDA_ERROR(cudaMemset(this->param_.d_data_, val, sizeof(T) * this->param_.number_of_elements_), "cudaMemset");
      }
      else
      {
        CHECK_CUDA_ERROR(micm::cuda::FillCudaMatrix<T>(this->param_, val), "FillCudaMatrix");
      }
    }

    CudaMatrixParam AsDeviceParam() const
    {
      return this->param_;
    }
  };
}  // namespace micm
