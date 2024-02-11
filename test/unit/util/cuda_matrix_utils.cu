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

      int SquareDriver(CudaVectorMatrixParam& param)
      {
          Square<<<param.num_elements_, 1>>>(param.d_data_, param.num_elements_);
      }
    }
}
