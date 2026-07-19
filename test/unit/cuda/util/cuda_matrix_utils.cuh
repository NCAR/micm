#include <micm/cuda/util/cuda_param.hpp>
#include <micm/util/types.hpp>

namespace micm::cuda
{
  void SquareDriver(CudaMatrixParam& param);
  void AddOneDriver(CudaMatrixParam& param);
  void SparseMatrixAddOneElementDriver(
      CudaMatrixParam& param,
      Index elem_id,
      Index grid_id,
      const Index cuda_matrix_vector_length);
  void DenseMatrixAddOneElementDriver(
      CudaMatrixParam& param,
      Index row_id,
      Index col_id,
      const Index cuda_matrix_vector_length);
}  // namespace micm::cuda
