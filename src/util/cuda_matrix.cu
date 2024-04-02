#include <cuda_runtime.h>

#include <micm/util/cuda_matrix.cuh>
#include <vector>

namespace micm
{
  namespace cuda
  {
    int MallocVector(CudaMatrixParam& param, std::size_t number_of_elements)
    {
      param.number_of_elements_ = number_of_elements;
      return cudaMalloc(&(param.d_data_), sizeof(double) * number_of_elements);
    }
    int FreeVector(CudaMatrixParam& param)
    {
      return cudaFree(param.d_data_);
    }
    int CopyToDevice(CudaMatrixParam& param, std::vector<double>& h_data)
    {
      return cudaMemcpy(param.d_data_, h_data.data(), sizeof(double) * param.number_of_elements_, cudaMemcpyHostToDevice);
    }
    int CopyToHost(CudaMatrixParam& param, std::vector<double>& h_data)
    {
      cudaDeviceSynchronize();
      return cudaMemcpy(h_data.data(), param.d_data_, sizeof(double) * param.number_of_elements_, cudaMemcpyDeviceToHost);
    }
    int CopyToDeviceFromDevice(CudaMatrixParam& vectorMatrixDest, const CudaMatrixParam& vectorMatrixSrc)
    {
      return cudaMemcpy(
          vectorMatrixDest.d_data_,
          vectorMatrixSrc.d_data_,
          sizeof(double) * vectorMatrixSrc.number_of_elements_,
          cudaMemcpyDeviceToDevice);
    }
  }  // namespace cuda
}  // namespace micm
