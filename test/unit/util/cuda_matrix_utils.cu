#include "cuda_matrix_utils.cuh"

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

      void SquareDriver(CudaVectorMatrixParam& param)
      {
          Square<<<param.number_of_elements_, 1>>>(param.d_data_, param.number_of_elements_);
      }

      __global__ void AddOne(double* d_data, std::size_t num_elements)
      {
          std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid < num_elements)
          {
              d_data[tid] += 1.0;
          }
      }

      void AddOneDriver(CudaVectorMatrixParam& param)
      {
        AddOne<<<param.number_of_elements_, BLOCK_SIZE>>>(param.d_data_, param.number_of_elements_);
      }
    }
}
