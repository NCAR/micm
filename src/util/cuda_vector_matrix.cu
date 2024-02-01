#include <vector>
//#include <micm/util/cuda_vector_matrix.cuh>

namespace micm
{
  namespace cuda
  {
    int MallocVector(double *d_data, std::size_t num_elements)
    {
      return cudaMalloc(&d_data, sizeof(double) * num_elements);
    }
    int FreeVector(double *d_data)
    {
      return cudaFree(d_data);
    }
    int CopyToDevice(double *d_data, const double *h_data, std::size_t num_elements)
    {
      return cudaMemcpy(d_data, h_data, sizeof(double) * num_elements, cudaMemcpyHostToDevice);
    }
    int CopyToHost(double* d_data, double *h_data, std::size_t num_elements)
    {
      return cudaMemcpy(h_data, d_data, sizeof(double) * num_elements, cudaMemcpyDeviceToHost);
    }
  }
}
