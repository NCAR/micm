#include <micm/util/cuda_param.hpp>
#include <micm/util/cuda_vector_matrix.cuh>
#include <micm/util/sparse_matrix.hpp>
#include <type_traits>

namespace micm
{
  template<class T, class OrderingPolicy>
  class CudaSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   private:
    CudaVectorMatrixParam param_;

   public:
    CudaSparseMatrix() = default;

    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(std::is_same_v<T, double>)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
      micm::cuda::MallocVector(param_, this->data_.size());
    }
    CudaSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : SparseMatrix<T, OrderingPolicy>(builder)
    {
    }

    CudaSparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder) requires(
        std::is_same_v<T, double>)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      micm::cuda::MallocVector(param_, this->data_.size());
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
      micm::cuda::MallocVector(param_, this->data_.size());
      micm::cuda::CopyToDeviceFromDevice(param_, other.param_);
    }

    CudaSparseMatrix(const CudaSparseMatrix& other)
        : SparseMatrix<T, OrderingPolicy>(other)
    {
    }

    CudaSparseMatrix(CudaSparseMatrix&& other) noexcept
        : SparseMatrix<T, OrderingPolicy>(SparseMatrix<T, OrderingPolicy>::create(other.number_of_blocks_))
    {
      this->data_ = std::move(other.data_);
      this->param_ = std::move(other.param_);
    }

    CudaSparseMatrix& operator=(const CudaSparseMatrix& other)
    {
      return *this = CudaSparseMatrix(other);
    }

    CudaSparseMatrix& operator=(CudaSparseMatrix&& other) noexcept
    {
      std::swap(this->data_, other.data_);
      std::swap(this->param_, other.param_);
      std::swap(this->number_of_blocks_, other.number_of_blocks_);
      std::swap(this->row_ids_, other.row_ids_);
      std::swap(this->row_start_, other.row_start_);
    }

    ~CudaSparseMatrix() requires(std::is_same_v<T, double>)
    {
      micm::cuda::FreeVector(param_);
    }

    ~CudaSparseMatrix()
    {
    }

    int CopyToDevice()
    {
      return micm::cuda::CopyToDevice(param_, this->data_);
    }
    int CopyToHost()
    {
      return micm::cuda::CopyToHost(param_, this->data_);
    }
    CudaVectorMatrixParam AsDeviceParam()
    {
      return this->param_;
    }
  };
}  // namespace micm
