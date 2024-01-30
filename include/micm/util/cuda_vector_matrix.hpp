#include <type_traits>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/cuda_vector_matrix.cuh>

namespace micm {

    template<class T, std::size_t L = DEFAULT_VECTOR_SIZE>
    class CudaVectorMatrix : public VectorMatrix<T, L> {
    private:
        double* d_data_;

    public:
        CudaVectorMatrix() requires(std::is_same_v<T, double>)
         : VectorMatrix<T,L>()
        {
          micm::cuda::malloc_vector(d_data_, static_cast<unsigned int>(0));
        }
        CudaVectorMatrix()
         : VectorMatrix<T,L>()
        {}

        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(x_dim, y_dim)
        {
          micm::cuda::malloc_vector(d_data_, this->data_.size());
        }
        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim)
          : VectorMatrix<T,L>(x_dim, y_dim)
        {}

        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(x_dim, y_dim, initial_value)
        {
          micm::cuda::malloc_vector(d_data_, this->data_.size());
        }
        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
          : VectorMatrix<T,L>(x_dim, y_dim, initial_value)
        {}

        CudaVectorMatrix(const std::vector<std::vector<T>> other) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(other)
        {
          micm::cuda::malloc_vector(d_data_, this->data_.size());
        }
        CudaVectorMatrix(const std::vector<std::vector<T>> other)
          : VectorMatrix<T,L>(other)
        {}

        ~CudaVectorMatrix() requires(std::is_same_v<T, double>)
        {
          micm::cuda::free_vector(d_data_);
        }

        int CopyToDevice()
        {
          static_assert(std::is_same_v<T, double>);
          return micm::cuda::copy_to_device(d_data_, this->data_.data(), this->AsVector().size());
        }
        int GetFromDevice()
        {
          static_assert(std::is_same_v<T, double>);
          return micm::cuda::copy_to_host(d_data_, this->data_.data(), this->AsVector().size());
        }
    };
}
