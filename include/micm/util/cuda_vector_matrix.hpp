#include <type_traits>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/cuda_vector_matrix.cuh>

namespace micm {

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
    * CUDA functionality requires T to be of type double, otherwise this
    * behaves similarily to VectorMatrix.
    */
    template<class T, std::size_t L = DEFAULT_VECTOR_SIZE>
    class CudaVectorMatrix : public VectorMatrix<T, L> {
    private:
        /// @brief The device pointer (handle) to the allocated memory on the target device.
        CudaVectorMatrixParam vector_matrix_param_;

    public:
        CudaVectorMatrix() requires(std::is_same_v<T, double>)
         : VectorMatrix<T,L>()
        {
          micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
        }
        CudaVectorMatrix()
         : VectorMatrix<T,L>()
        {}

        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(x_dim, y_dim)
        {
          micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
        }
        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim)
          : VectorMatrix<T,L>(x_dim, y_dim)
        {}

        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(x_dim, y_dim, initial_value)
        {
          micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
        }
        CudaVectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
          : VectorMatrix<T,L>(x_dim, y_dim, initial_value)
        {}

        CudaVectorMatrix(const std::vector<std::vector<T>> other) requires(std::is_same_v<T, double>)
          : VectorMatrix<T,L>(other)
        {
          micm::cuda::MallocVector(vector_matrix_param_, this->data_.size());
        }
        CudaVectorMatrix(const std::vector<std::vector<T>> other)
          : VectorMatrix<T,L>(other)
        {}

        ~CudaVectorMatrix() requires(std::is_same_v<T, double>)
        {
          micm::cuda::FreeVector(vector_matrix_param_);
        }

        int CopyToDevice()
        {
          static_assert(std::is_same_v<T, double>);
          return micm::cuda::CopyToDevice(vector_matrix_param_, this->AsVector());
        }
        int CopyToHost()
        {
          static_assert(std::is_same_v<T, double>);
          return micm::cuda::CopyToHost(vector_matrix_param_, this->AsVector());
        }
        CudaVectorMatrixParam AsDeviceParam()
        {
          return CudaVectorMatrixParam{ vector_matrix_param_.d_data_, vector_matrix_param_.num_elements_ };
        }
    };
}
