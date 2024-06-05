/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once
 
#include <cublas_v2.h>
#include <cuda_runtime.h> 
#include <string>

#define CHECK_CUDA_ERROR(err, msg)   micm::cuda::CheckCudaError(err, __FILE__, __LINE__, msg)
#define CHECK_CUBLAS_ERROR(err, msg) micm::cuda::CheckCublasError(err, __FILE__, __LINE__, msg)

namespace micm
{
  namespace cuda
  {
    /// @brief Checks for CUDA errors and prints error message if any
    /// @param err Error code to check
    /// @param file File where error occurred
    /// @param line Line number where error occurred
    /// @param str Additional string to print with error message
    void CheckCudaError(cudaError_t err, const char* file, int line, std::string str);

    /// @brief Checks for cuBLAS errors and prints error message if any
    /// @param err Error code to check
    /// @param file File where error occurred
    /// @param line Line number where error occurred
    /// @param str Additional string to print with error message
    void CheckCublasError(cublasStatus_t err, const char* file, int line, std::string str);
  }  // namespace cuda
}  // namespace micm