// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <map>
#include <memory>
#include <mutex>
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

    /// @brief Get the cuBLAS handle for the current device
    cublasHandle_t& GetCublasHandle();

    /// @brief Define a functor for the cudaStream unique pointer deleter
    struct CudaStreamDeleter
    {
      void operator()(cudaStream_t* cuda_stream) const
      {
        if (cuda_stream != nullptr)
        {
          CHECK_CUDA_ERROR(cudaStreamDestroy(*cuda_stream), "CUDA stream finalization failed");
          delete cuda_stream;
          cuda_stream = nullptr;
        }
      }
    };

    /// @brief Define the smart pointer type using the functor for the custom deleter
    using CudaStreamPtr = std::unique_ptr<cudaStream_t, CudaStreamDeleter>;

    /// @brief Define a functor for the cudaEvent unique pointer deleter
    struct CudaEventDeleter
    {
      void operator()(cudaEvent_t* cuda_event) const
      {
        if (cuda_event != nullptr)
        {
          CHECK_CUDA_ERROR(cudaEventDestroy(*cuda_event), "CUDA event finalization failed");
          delete cuda_event;
          cuda_event = nullptr;
        }
      }
    };

    /// @brief Define the smart pointer type using the functor for the custom deleter
    using CudaEventPtr = std::unique_ptr<cudaEvent_t, CudaEventDeleter>;

    /// @brief Singleton class to manage CUDA streams
    class CudaStreamSingleton
    {
     public:
      ~CudaStreamSingleton() = default;

      CudaStreamSingleton(const CudaStreamSingleton&) = delete;

      CudaStreamSingleton& operator=(const CudaStreamSingleton&) = delete;

      // Get the only one instance of the singleton class
      static CudaStreamSingleton& GetInstance();

      // Get the CUDA stream given a stream ID
      cudaStream_t& GetCudaStream(std::size_t stream_id);

      // Get the CUDA event given an event ID
      cudaEvent_t& GetCudaEvent(std::size_t event_id);

      // Empty the map variable to clean up all CUDA streams
      void CleanUp();

     private:
      // Private constructor to prevent direct instantiation
      CudaStreamSingleton() = default;

      // Create a CUDA stream and return a unique pointer to it
      CudaStreamPtr CreateCudaStream();

      // Create a CUDA event and return a unique pointer to it
      CudaEventPtr CreateCudaEvent();

      std::map<int, CudaStreamPtr> cuda_streams_map_;

      std::map<int, CudaEventPtr> cuda_events_map_;
      
      std::mutex mutex_;
    };
  }  // namespace cuda
}  // namespace micm