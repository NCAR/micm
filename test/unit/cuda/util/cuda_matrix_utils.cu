#include "cuda_matrix_utils.cuh"

#include <micm/util/types.hpp>

#include <cstdio>
#include <iostream>

namespace micm::cuda
{
  __global__ void Square(Real* d_data, Index num_elements)
  {
    Index tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < num_elements)
    {
      d_data[tid] *= d_data[tid];
    }
  }

  void SquareDriver(CudaMatrixParam& param)
  {
    Square<<<param.number_of_elements_, 1>>>(param.d_data_, param.number_of_elements_);
  }

  __global__ void AddOne(Real* d_data, Index num_elements)
  {
    Index tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < num_elements)
    {
      d_data[tid] += 1;
    }
  }

  __global__ void SparseMatrixAddOneElement(
      Real* d_data,
      const Index number_of_grid_cells,
      const Index number_of_non_zeros,
      Index elem_id,
      Index grid_id,
      const Index cuda_matrix_vector_length)
  {
    Index tid = blockIdx.x * blockDim.x + threadIdx.x;
    Index local_tid = tid % cuda_matrix_vector_length;
    Index group_id = std::floor(static_cast<Real>(tid) / cuda_matrix_vector_length);
    const Index offset =
        group_id * number_of_non_zeros * cuda_matrix_vector_length + elem_id * cuda_matrix_vector_length + local_tid;
    if (tid == grid_id)
    {
      d_data[offset] += 1;
    }
  }

  __global__ void DenseMatrixAddOneElement(
      Real* d_data,
      const Index number_of_columns,
      Index row_id,
      Index col_id,
      const Index cuda_matrix_vector_length)
  {
    Index tid = blockIdx.x * blockDim.x + threadIdx.x;
    Index local_tid = tid % cuda_matrix_vector_length;
    Index group_id = std::floor(static_cast<Real>(tid) / cuda_matrix_vector_length);
    const Index offset =
        group_id * number_of_columns * cuda_matrix_vector_length + col_id * cuda_matrix_vector_length + local_tid;
    if (tid == row_id)
    {
      d_data[offset] += 1;
    }
  }

  void AddOneDriver(CudaMatrixParam& param)
  {
    AddOne<<<param.number_of_elements_, BLOCK_SIZE>>>(param.d_data_, param.number_of_elements_);
  }

  void SparseMatrixAddOneElementDriver(
      CudaMatrixParam& param,
      Index elem_id,
      Index grid_id,
      const Index cuda_matrix_vector_length)
  {
    if (grid_id >= param.number_of_grid_cells_)
    {
      throw std::runtime_error("grid_id out of bounds in SparseMatrixAddOneElementDriver");
    }
    const Index number_of_groups =
        std::ceil(static_cast<Real>(param.number_of_grid_cells_) / cuda_matrix_vector_length);
    const Index number_of_non_zeros = param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;
    const Index number_of_blocks = (param.number_of_grid_cells_ + 31) / 32;
    SparseMatrixAddOneElement<<<number_of_blocks, 32>>>(
        param.d_data_, param.number_of_grid_cells_, number_of_non_zeros, elem_id, grid_id, cuda_matrix_vector_length);
  }

  void DenseMatrixAddOneElementDriver(
      CudaMatrixParam& param,
      Index row_id,
      Index col_id,
      const Index cuda_matrix_vector_length)
  {
    if (row_id >= param.number_of_grid_cells_)
    {
      throw std::runtime_error("row_id out of bounds in DenseMatrixAddOneElementDriver");
    }
    const Index number_of_groups =
        std::ceil(static_cast<Real>(param.number_of_grid_cells_) / cuda_matrix_vector_length);
    const Index number_of_columns = param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;
    const Index number_of_blocks = (param.number_of_grid_cells_ + 31) / 32;
    DenseMatrixAddOneElement<<<number_of_blocks, 32>>>(
        param.d_data_, number_of_columns, row_id, col_id, cuda_matrix_vector_length);
  }
}  // namespace micm::cuda
