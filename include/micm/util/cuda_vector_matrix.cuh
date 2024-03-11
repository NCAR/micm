#include <vector>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// @brief Allocate memory on device
    /// @param vectorMatrix Reference to struct containing information about allocated memory
    /// @param num_elements Requested number of elements to allocate
    /// @returns Error code from allocating data on the device, if any
    int MallocVector(CudaVectorMatrixParam& vectorMatrix, std::size_t num_elements);

    /// @brief Free memory allocated on device
    /// @param vectorMatrix Struct containing allocated device memory
    /// @returns Error code from free-ing data on device, if any
    int FreeVector(CudaVectorMatrixParam& vectorMatrix);

    /// @brief Copies data from the host to the device
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param h_data Host data to copy from
    /// @returns Error code from copying to device from the host, if any
    int CopyToDevice(CudaVectorMatrixParam& vectorMatrix, std::vector<double>& h_data);

    /// @brief Copies data from the device to the host
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param h_data Host data to copy data to
    /// @returns Error code from copying from the device to the host, if any
    int CopyToHost(CudaVectorMatrixParam& vectorMatrix, std::vector<double>& h_data);

    /// @brief Copies data to the destination device memory block from the source device memory block
    /// @param vectorMatrixDest Struct containing allocated destination device memory to copy to
    /// @param vectorMatrixSrc Struct containing allocated source device memory to copy from
    /// @returns Error code from copying to destination device memory from source device memory, if any
    int CopyToDeviceFromDevice(CudaVectorMatrixParam& vectorMatrixDest, const CudaVectorMatrixParam& vectorMatrixSrc);
  }
}
