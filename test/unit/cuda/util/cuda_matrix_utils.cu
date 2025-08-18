#include "cuda_matrix_utils.cuh"
#include <cstdio>
#include <iostream>

namespace micm
{
  namespace cuda
  {
    __global__ void Square(double* d_data, std::size_t num_elements)
    {
      std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      if (tid < num_elements)
      {
        d_data[tid] *= d_data[tid];
      }
    }

    void SquareDriver(CudaMatrixParam& param)
    {
      Square<<<param.number_of_elements_, 1>>>(param.d_data_, param.number_of_elements_);
    }

    __global__ void AddOne(double* d_data, std::size_t num_elements)
    {
      std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      if (tid < num_elements)
      {
        d_data[tid] += 1;
      }
    }

    __global__ void AddOneElement(double* d_data, const std::size_t number_of_non_zeros, std::size_t elem_id, std::size_t grid_id, const std::size_t cuda_matrix_vector_length)
    {
      std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
      std::size_t local_tid = tid % cuda_matrix_vector_length;
      std::size_t group_id = std::floor(static_cast<double>(tid) / cuda_matrix_vector_length);
      const std::size_t offset = group_id * number_of_non_zeros * cuda_matrix_vector_length + elem_id * cuda_matrix_vector_length + local_tid;
      if (tid == grid_id)
      {
        d_data[offset] += 1;
      }
    }

    void AddOneDriver(CudaMatrixParam& param)
    {
      AddOne<<<param.number_of_elements_, BLOCK_SIZE>>>(param.d_data_, param.number_of_elements_);
    }

    void AddOneElementDriver(CudaMatrixParam& param, std::size_t elem_id, std::size_t grid_id, const std::size_t cuda_matrix_vector_length)
    {
      const std::size_t number_of_non_zeros = param.number_of_elements_ / param.number_of_grid_cells_; 
      AddOneElement<<<param.number_of_grid_cells_, 32>>>(param.d_data_, number_of_non_zeros, elem_id, grid_id, cuda_matrix_vector_length);
    }

  }  // namespace cuda
}  // namespace micm
