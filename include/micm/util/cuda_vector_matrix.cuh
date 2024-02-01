#include <vector>

namespace micm
{
  namespace cuda
  {
    int MallocVector(double *d_data, std::size_t num_elements);

    int FreeVector(double *d_data);

    int CopyToDevice(double *d_data, const double *h_data, std::size_t num_elements);

    int CopyToHost(double *d_data, double *h_data, std::size_t num_elements);
  }
}
