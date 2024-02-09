#include <vector>
#include <micm/util/cuda_vector_matrix.cuh>

namespace micm
{
  namespace cuda
  {
    int MallocVector(CudaVectorMatrixParam& param, std::size_t num_elements)
    {
      param.num_elements_ = num_elements;
      return cudaMalloc(&(param.d_data_), sizeof(double) * num_elements);
    }
    int FreeVector(CudaVectorMatrixParam& param)
    {
      return cudaFree(param.d_data_);
    }
    int CopyToDevice(CudaVectorMatrixParam& param, std::vector<double>& h_data)
    {
      return cudaMemcpy(param.d_data_, h_data.data(), sizeof(double) * param.num_elements_, cudaMemcpyHostToDevice);
    }
    int CopyToHost(CudaVectorMatrixParam& param, std::vector<double>& h_data)
    {
      return cudaMemcpy(h_data.data(), param.d_data_, sizeof(double) * param.num_elements_, cudaMemcpyDeviceToHost);
    }
  }
}
