// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

#include <cuda_runtime.h>

#include <vector>

namespace micm
{
  namespace cuda
  {
    /// @brief Allocate memory on device
    /// @param vectorMatrix Reference to struct containing information about allocated memory
    /// @param num_elements Requested number of elements to allocate
    /// @returns Error code from allocating data on the device, if any
    template<typename T>
    cudaError_t MallocVector(CudaMatrixParam& vectorMatrix, std::size_t num_elements);

    /// @brief Free memory allocated on device
    /// @param vectorMatrix Struct containing allocated device memory
    /// @returns Error code from free-ing data on device, if any
    cudaError_t FreeVector(CudaMatrixParam& vectorMatrix);

    /// @brief Sets all elements in a CUDA matrix to their current value or a specified value, whichever is greater
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param val Value to compare with each element in the matrix
    template<typename T>
    cudaError_t MatrixMax(CudaMatrixParam& vectorMatrix, T val);

    /// @brief Sets all elements in a CUDA matrix to their current value or a specified value, whichever is lesser
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param val Value to compare with each element in the matrix
    template<typename T>
    cudaError_t MatrixMin(CudaMatrixParam& vectorMatrix, T val);

    /// @brief Copies data from the host to the device
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param h_data Host data to copy from
    /// @returns Error code from copying to device from the host, if any
    template<typename T>
    cudaError_t CopyToDevice(CudaMatrixParam& vectorMatrix, const std::vector<T>& h_data);

    /// @brief Copies data from the device to the host
    /// @param vectorMatrix Struct containing allocated device memory
    /// @param h_data Host data to copy data to
    /// @returns Error code from copying from the device to the host, if any
    template<typename T>
    cudaError_t CopyToHost(CudaMatrixParam& vectorMatrix, std::vector<T>& h_data);

    /// @brief Copies data to the destination device memory block from the source device memory block
    /// @param vectorMatrixDest Struct containing allocated destination device memory to copy to
    /// @param vectorMatrixSrc Struct containing allocated source device memory to copy from
    /// @returns Error code from copying to destination device memory from source device memory, if any
    template<typename T>
    cudaError_t CopyToDeviceFromDevice(CudaMatrixParam& vectorMatrixDest, const CudaMatrixParam& vectorMatrixSrc);

    /// @brief Fills a CUDA matrix with a specified value
    /// @param param Struct containing allocated device memory
    /// @param val Value to fill the matrix with
    template<typename T>
    cudaError_t FillCudaMatrix(CudaMatrixParam& param, T val);
  }  // namespace cuda
}  // namespace micm
