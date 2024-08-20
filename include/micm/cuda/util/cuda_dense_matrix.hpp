// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/util/error.hpp>
#include <micm/util/vector_matrix.hpp>

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <type_traits>

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
   * Copy function only copies the device data from one CUDA dense matrix
   * to the other (no change of the host data), assuming that the device memory
   * has been allocated correctly. A check is done before doing the copy
   * to make sure that both matrices have the same size.
   *
   * CUDA functionality requires T to be of type double, otherwise this
   * behaves similarily to VectorMatrix.
   */

  /// Concept for Cuda Matrix
  template<typename MatrixType>
  concept CudaMatrix = requires(MatrixType t)
  {
    {
      t.CopyToDevice()
      } -> std::same_as<void>;
    {
      t.CopyToHost()
      } -> std::same_as<void>;
    {
      t.AsDeviceParam()
      } -> std::same_as<CudaMatrixParam>;
  };

  template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class CudaDenseMatrix : public VectorMatrix<T, L>
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = CudaDenseMatrix<int, L>;
    using value_type = T;

   private:
    /// @brief The device pointer (handle) to the allocated memory on the target device.
    CudaMatrixParam param_;

   public:
    CudaDenseMatrix() requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>()
    {
      this->param_.number_of_grid_cells_ = 0;
      this->param_.number_of_elements_ = this->data_.size();
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
    }
    CudaDenseMatrix()
        : VectorMatrix<T, L>()
    {
    }

    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
      this->param_.number_of_elements_ = this->data_.size();
      this->param_.number_of_grid_cells_ = x_dim;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
    }
    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
    }

    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
      this->param_.number_of_elements_ = this->data_.size();
      this->param_.number_of_grid_cells_ = x_dim;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice(this->param_, this->data_), "cudaMemcpyHostToDevice");
    }
    CudaDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
    }

    CudaDenseMatrix(const std::vector<std::vector<T>> other) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(other)
    {
      this->param_.number_of_grid_cells_ = 0;
      this->param_.number_of_elements_ = 0;
      for (const auto& inner_vector : other)
      {
        this->param_.number_of_elements_ += inner_vector.size();
      }
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
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
      this->param_.number_of_elements_ = other.param_.number_of_elements_;
      this->param_.number_of_grid_cells_ = other.param_.number_of_grid_cells_;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
    }

    CudaDenseMatrix(const CudaDenseMatrix& other)
        : VectorMatrix<T, L>(other)
    {
    }

    CudaDenseMatrix(CudaDenseMatrix&& other) noexcept
        : VectorMatrix<T, L>(other)
    {
      this->param_.d_data_ = nullptr;
      std::swap(this->param_, other.param_);
    }

    CudaDenseMatrix& operator=(const CudaDenseMatrix& other)
    {
      VectorMatrix<T, L>::operator=(other);
      if (this->param_.d_data_ != nullptr)
        CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), "cudaFree");
      this->param_ = other.param_;
      CHECK_CUDA_ERROR(micm::cuda::MallocVector(this->param_, this->param_.number_of_elements_), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
      return *this;
    }

    CudaDenseMatrix& operator=(CudaDenseMatrix&& other) noexcept
    {
      if (this != &other)
      {
        VectorMatrix<T, L>::operator=(other);
        std::swap(this->param_, other.param_);
      }
      return *this;
    }

    ~CudaDenseMatrix()
    {
    }

    ~CudaDenseMatrix() requires(std::is_same_v<T, double>)
    {
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(this->param_), "cudaFree");
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

    /// @brief For each element in the VectorMatrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant.
    /// @param alpha The scaling scalar to apply to the VectorMatrix x
    /// @param x The input VectorMatrix
    /// @return 0 if successful, otherwise an error code
    void Axpy(const double alpha, const CudaDenseMatrix<T, L>& x)
    {
      const int incx = 1;  // increment for the elements of x
      const int incy = 1;  // increment for the elements of y
      static_assert(std::is_same_v<T, double>);
      CHECK_CUBLAS_ERROR(
          cublasDaxpy(
              micm::cuda::GetCublasHandle(),
              x.param_.number_of_elements_,
              &alpha,
              x.param_.d_data_,
              incx,
              this->param_.d_data_,
              incy),
          "CUBLAS Daxpy operation failed...");
    }

    // Copy the device data from the other Cuda dense matrix into this one
    void Copy(const CudaDenseMatrix& other)
    {
      if (other.param_.number_of_elements_ != this->param_.number_of_elements_)
      {
        throw std::runtime_error("Both CUDA dense matrices must have the same size.");
      }
      CHECK_CUDA_ERROR(micm::cuda::CopyToDeviceFromDevice(this->param_, other.param_), "cudaMemcpyDeviceToDevice");
    }

    /// @brief Set every matrix element to a given value on the GPU
    /// @param val Value to set each element to
    void Fill(T val)
    {
      if constexpr (std::is_same_v<T, int>)
      {
        // the cudaMemset function only works for integer types and is an asynchronous function: https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1gf7338650f7683c51ee26aadc6973c63a
        CHECK_CUDA_ERROR(cudaMemset(this->param_.d_data_, val, sizeof(double) * this->param_.number_of_elements_), "cudaMemset");
      }
      else
      {
        CHECK_CUDA_ERROR(micm::cuda::FillCudaMatrix<T>(this->param_, val), "FillCudaMatrix");
      }
    }

  };  // class CudaDenseMatrix
}  // namespace micm
