#include <micm/util/cuda_vector_matrix.cuh>
#include <micm/util/vector_matrix.hpp>
#include <type_traits>

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
  class CudaVectorMatrix : public VectorMatrix<T, L>
  {
   private:
    /// @brief The device pointer (handle) to the allocated memory on the target device.
    CudaVectorMatrixParam vector_matrix_param_;
    /// @brief The handle to the CUBLAS library
    cublasHandle_t handle_ = NULL;

   public:
    CudaVectorMatrix() requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>()
    {
      micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
    }
    CudaVectorMatrix()
        : VectorMatrix<T, L>()
    {
    }

    CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
      micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
      this->vector_matrix_param_.number_of_grid_cells_ = x_dim;
    }
    CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
    }

    CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
      micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
      this->vector_matrix_param_.number_of_grid_cells_ = x_dim;
      if (this->handle_ == NULL)
      {
        cublasStatus_t stat = cublasCreate(&(this->handle_));
        if (stat != CUBLAS_STATUS_SUCCESS)
        {
          std::cout << stat << std::endl;
          throw std::runtime_error("CUBLAS initialization failed.");
        }
      }
    }
    CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
    }

    CudaVectorMatrix(const std::vector<std::vector<T>> other) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(other)
    {
      micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
    }

    CudaVectorMatrix(const std::vector<std::vector<T>> other)
        : VectorMatrix<T, L>(other)
    {
    }

    CudaVectorMatrix(const CudaVectorMatrix& other) requires(std::is_same_v<T, double>)
        : VectorMatrix<T, L>(other.x_dim_, other.y_dim_)
    {
      this->data_ = other.data_;
      micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
      micm::cuda::CopyToDeviceFromDevice(vector_matrix_param_, other.vector_matrix_param_);
    }

    CudaVectorMatrix(const CudaVectorMatrix& other)
        : VectorMatrix<T, L>(other.x_dim_, other.y_dim_)
    {
      this->data_ = other.data_;
    }

    CudaVectorMatrix(CudaVectorMatrix&& other) noexcept
        : VectorMatrix<T, L>(other.x_dim_, other.y_dim_)
    {
      this->data_ = std::move(other.data_);
      this->vector_matrix_param_ = std::move(other.vector_matrix_param_);
    }

    CudaVectorMatrix& operator=(const CudaVectorMatrix& other)
    {
      return *this = CudaVectorMatrix(other);
    }

    CudaVectorMatrix& operator=(CudaVectorMatrix&& other) noexcept
    {
      std::swap(this->data_, other.data_);
      std::swap(this->vector_matrix_param_, other.vector_matrix_param_);
      this->x_dim_ = other.x_dim_;
      this->y_dim_ = other.y_dim_;
      return *this;
    }

    ~CudaVectorMatrix()
    {
    }

    ~CudaVectorMatrix() requires(std::is_same_v<T, double>)
    {
      micm::cuda::FreeVector(this->vector_matrix_param_);
      if (this->handle_ != NULL)
        cublasDestroy(this->handle_);
    }

    int CopyToDevice()
    {
      static_assert(std::is_same_v<T, double>);
      return micm::cuda::CopyToDevice(vector_matrix_param_, this->data_);
    }
    int CopyToHost()
    {
      static_assert(std::is_same_v<T, double>);
      return micm::cuda::CopyToHost(vector_matrix_param_, this->data_);
    }
    CudaVectorMatrixParam AsDeviceParam() const
    {
      return this->vector_matrix_param_;
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
    void Axpy(const double alpha, const CudaVectorMatrix<T, L>& x, const int incx, const int incy)
    {
      static_assert(std::is_same_v<T, double>);
      cublasStatus_t stat = cublasDaxpy(
          this->handle_,
          x.vector_matrix_param_.number_of_elements_,
          &alpha,
          x.vector_matrix_param_.d_data_,
          incx,
          this->vector_matrix_param_.d_data_,
          incy);
      if (stat != CUBLAS_STATUS_SUCCESS)
      {
        throw std::runtime_error("CUBLAS Daxpy operation failed.");
      }
    }
  };
}  // namespace micm
