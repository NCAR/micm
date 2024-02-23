#include <micm/util/cuda_vector_matrix.cuh>
#include <vector>

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
    int CopyToDeviceFromDevice(CudaVectorMatrixParam& vectorMatrixDest, const CudaVectorMatrixParam& vectorMatrixSrc)
    {
      return cudaMemcpy(
          vectorMatrixDest.d_data_,
          vectorMatrixSrc.d_data_,
          sizeof(double) * vectorMatrixSrc.num_elements_,
          cudaMemcpyDeviceToDevice);
    }
  }  // namespace cuda
}  // namespace micm
