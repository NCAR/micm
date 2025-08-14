#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    void SquareDriver(CudaMatrixParam& param);
    void AddOneDriver(CudaMatrixParam& param);
    void AddOneElementDriver(CudaMatrixParam& param, std::size_t elem_id, std::size_t grid_id, const std::size_t cuda_matrix_vector_length);
  }  // namespace cuda
}  // namespace micm
