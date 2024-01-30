#include <vector>

namespace micm
{
  namespace cuda
  {
    int malloc_vector(double *d_data, std::size_t num_elements);

    int free_vector(double *d_data);

    int copy_to_device(double *d_data, const double *h_data, std::size_t num_elements);

    int copy_to_host(double *d_data, double *h_data, std::size_t num_elements);
  }
}
