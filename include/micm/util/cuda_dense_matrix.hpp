// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/cuda_matrix.cuh>
#include <micm/util/vector_matrix.hpp>
#include <type_traits>
#include <cuda_runtime.h>
#include "cublas_v2.h"

namespace micm
{

  /**
   * @brief Provides a CUDA implemtation to the VectorMatrix functionality.
   *
   * This class provides the VectorMatrix API but allows for operations to
   * be performed via CUDA with the requirement that the caller explicity
   * move data to and from the device.
   *
   * After performing operations with ForEach, the caller must decide
   * when to syncronize
   * the host data with GetFromDevice() and any modification of host data
   * including initialization must be followed by CopyToDevice() otherwise
   * host and device data will be out of sync.
   *
   * Copy/Move constructors/assignment operators are non-synchronizing
   * operators/constructors so if device and host data is desynchronized,
   * the copies and moved matrices will remain desynchronized.
   *
   * CUDA functionality requires T to be of type double, otherwise this
   * behaves similarily to VectorMatrix.
   */
  template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class CudaDenseMatrix : public VectorMatrix<T, L>
  {
   private:
    /// @brief The device pointer (handle) to the allocated memory on the target device.
    CudaMatrixParam param_;
    /// @brief The handle to the CUBLAS library
    cublasHandle_t handle_ = NULL;

   public:
    CudaDenseMatrix() requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>()
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
    }
    CudaDenseMatrix()
        : VectorMatrix<T, L>()
    {
    }

    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      this->param_.number_of_grid_cells_ = x_dim;
    }
    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
    }

    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      this->param_.number_of_grid_cells_ = x_dim;
      micm::cuda::CHECK_CUBLAS_ERROR(cublasCreate(&(this->handle_)), __FILE__, __LINE__, "CUBLAS initialization failed...");
    }
    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
    }

    CudaDenseMatrix(const std::vector<std::vector<T>> other) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(other)
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
    }

    CudaDenseMatrix(const std::vector<std::vector<T>> other)
        : VectorMatrix<T, L>(other)
    {
    }

    CudaDenseMatrix(const CudaDenseMatrix& other) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(other)
    {
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), __FILE__, __LINE__, "cudaMemcpyDeviceToDevice");
      micm::cuda::CHECK_CUBLAS_ERROR(cublasCreate(&(this->handle_)), __FILE__, __LINE__, "CUBLAS initialization failed...");
    }

    CudaDenseMatrix(const CudaDenseMatrix& other)
        : VectorMatrix<T, L>(other)
    {
    }

    CudaDenseMatrix(CudaDenseMatrix&& other) noexcept
        : VectorMatrix<T, L>(other)
    {
      std::swap(this->param_, other.param_);
      std::swap(this->handle_, other.handle_);
    }

    CudaDenseMatrix& operator=(const CudaDenseMatrix& other)
    {
      VectorMatrix<T, L>::operator=(other);
      this->param_ = other.param_;
      this->param_.d_data_ = nullptr;
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->data_.size()), __FILE__, __LINE__, "cudaMalloc");
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), __FILE__, __LINE__, "cudaMemcpyDeviceToDevice");
      micm::cuda::CHECK_CUBLAS_ERROR(cublasCreate(&(this->handle_)), __FILE__, __LINE__, "CUBLAS initialization failed...");
      return *this;
    }

    CudaDenseMatrix& operator=(CudaDenseMatrix&& other) noexcept
    {
      if (this != &other)
      {
        VectorMatrix<T, L>::operator=(other);
        std::swap(this->param_, other.param_);
        std::swap(this->handle_, other.handle_);
      }
      return *this;
    }

    ~CudaDenseMatrix()
    {
    }

    ~CudaDenseMatrix() requires(std::is_same_v<T, double>)
    {
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), __FILE__, __LINE__, "cudaFree");
      if (this->handle_ != NULL)
        cublasDestroy(this->handle_);
      this->param_.d_data_ = nullptr;
      this->handle_ = NULL;
    }

    void CopyToDevice()
    {
      static_assert(std::is_same_v<T, double>);
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToDevice(this->param_, this->data_), __FILE__, __LINE__, "cudaMemcpyHostToDevice");
    }

    void CopyToHost()
    {
      static_assert(std::is_same_v<T, double>);
      micm::cuda::CHECK_CUDA_ERROR(micm::cuda::CopyToHost(this->param_, this->data_), __FILE__, __LINE__, "cudaMemcpyDeviceToHost");
    }

    CudaMatrixParam AsDeviceParam() const
    {
      return this->param_;
    }

    cublasHandle_t AsCublasHandle() const
    {
      return this->handle_;
    }

    /// @brief For each element in the VectorMatrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant.
    /// @param alpha The scaling scalar to apply to the VectorMatrix x
    /// @param x The input VectorMatrix
    /// @param incx The increment for the elements of x
    /// @param incy The increment for the elements of y
    /// @return 0 if successful, otherwise an error code
    void Axpy(const double alpha, const CudaDenseMatrix<T, L>& x, const int incx, const int incy)
    {
      static_assert(std::is_same_v<T, double>);
      micm::cuda::CHECK_CUBLAS_ERROR(cublasDaxpy(this->handle_, x.param_.number_of_elements_, &alpha, x.param_.d_data_, incx, this->param_.d_data_, incy), __FILE__, __LINE__, "CUBLAS Daxpy operation failed...");
    }
  };
}  // namespace micm
