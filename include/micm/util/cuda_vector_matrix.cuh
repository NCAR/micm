#include <vector>

namespace micm
{
  namespace cuda
  {
    /// @brief Allocate memory on device
    /// @param d_data Pointer to desired allocated memory
    /// @param Requested number of elements to allocate
    /// @returns Error code from allocating data on the device, if any
    int MallocVector(double *d_data, std::size_t num_elements);

    /// @brief Free memory allocated on device
    /// @param Device pointer pointing to allocated data
    /// @returns Error code from free-ing data on device.
    int FreeVector(double *d_data);

    /// @brief Copies data from the host to the device
    /// @param d_data Device pointer to copy data to
    /// @param h_data Host data to copy from
    /// @param num_elements The number of elements to copy to the device
    /// @returns Error code from copying to device from the host, if any
    int CopyToDevice(double *d_data, const double *h_data, std::size_t num_elements);

    /// @brief Copies data from the device to the host
    /// @param d_data Device pointer to copy data from
    /// @param h_data Host data to copy data to
    /// @param num_elements The number of elements to copy from the device
    /// @returns Error code from copying from the device to the host, if any
    int CopyToHost(double *d_data, double *h_data, std::size_t num_elements);
  }
}
