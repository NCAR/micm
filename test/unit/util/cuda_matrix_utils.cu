#include "cuda_matrix_utils.cuh"

namespace micm
{
    namespace cuda
    {
      __global__ void MultiplyByIndex(double* d_data, std::size_t num_elements)
      {
          std::size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
          if (tid < num_elements)
          {
              d_data[tid] *= tid;
          }
      }

      int MutiplyByIndexDriver(CudaVectorMatrixParam& param)
      {
          MultiplyByIndex<<<param.num_elements_, 1>>>(param.d_data_, param.num_elements_);
      }
    }
}
