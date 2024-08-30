// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/util/internal_error.hpp>

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <map>
#include <memory>
#include <mutex>

namespace micm
{
  namespace cuda
  {
    void CheckCudaError(cudaError_t err, const char* file, int line, std::string str)
    {
      if (err != cudaSuccess)
      {
        std::string msg = std::string(cudaGetErrorString(err)) + " : " + str;
        ThrowInternalError(MicmInternalErrc::Cuda, file, line, msg.c_str());
      }
    }

    void CheckCublasError(cublasStatus_t err, const char* file, int line, std::string str)
    {
      if (err != CUBLAS_STATUS_SUCCESS)
      {
        std::string msg = std::to_string(err) + " : " + str;
        ThrowInternalError(MicmInternalErrc::Cublas, file, line, msg.c_str());
      }
    }

    /*
       The following functions are used to create and manage cublas handles
    */

    // Define a functor for the cublasHandle_t unique pointer deleter
    struct CublasHandleDeleter
    {
      void operator()(cublasHandle_t* handle) const
      {
        if (handle != nullptr)
        {
          CHECK_CUBLAS_ERROR(cublasDestroy(*handle), "CUBLAS finalization failed");
          delete handle;
        }
      }
    };

    // Define the smart pointer type using the functor for the custom deleter
    using CublasHandlePtr = std::unique_ptr<cublasHandle_t, CublasHandleDeleter>;

    // Create a cublas handle and return a unique pointer to it
    CublasHandlePtr CreateCublasHandle()
    {
      cublasHandle_t* handle = new cublasHandle_t;
      CHECK_CUBLAS_ERROR(cublasCreate(handle), "CUBLAS initialization failed...");
      return CublasHandlePtr(handle, CublasHandleDeleter());
    }

    // Get the cublas handle for the current device
    cublasHandle_t& GetCublasHandle()
    {
      static std::map<int, CublasHandlePtr> cublas_handles_map;
      static std::mutex mutex;
      int device_id;
      CHECK_CUDA_ERROR(cudaGetDevice(&device_id), "Failed to get device ID...");
      std::lock_guard<std::mutex> lock(mutex);  // lock the mutex and generate a new cublas handle below
      if (auto search = cublas_handles_map.find(device_id); search == cublas_handles_map.end())
      {
        cublas_handles_map[device_id] = std::move(CreateCublasHandle());
      }
      return *cublas_handles_map[device_id];
    }

    /*
        The following functions are used to create and manage cuda streams
    */

    // Define a functor for the cudaStream unique pointer deleter
    struct CudaStreamDeleter
    {
      void operator()(cudaStream_t* cuda_stream) const
      {
        if (cuda_stream != nullptr)
        {
          cudaStreamSynchronize(*cuda_stream);
          CHECK_CUDA_ERROR(cudaStreamDestroy(*cuda_stream), "CUDA stream finalization failed");
          delete cuda_stream;
        }
      }
    };

    // Define the smart pointer type using the functor for the custom deleter
    using CudaStreamPtr = std::unique_ptr<cudaStream_t, CudaStreamDeleter>;

    // Create a CUDA stream and return a unique pointer to it
    CudaStreamPtr CreateCudaStream()
    {
      cudaStream_t* cuda_stream = new cudaStream_t;
      CHECK_CUDA_ERROR(cudaStreamCreate(cuda_stream), "CUDA stream initialization failed...");
      return CudaStreamPtr(cuda_stream, CudaStreamDeleter());
    }

    // Get the CUDA stream given a stream ID
    cudaStream_t& GetCudaStream(std::size_t stream_id)
    {
      static std::map<int, CudaStreamPtr> cuda_streams_map;
      if (auto search = cuda_streams_map.find(stream_id); search == cuda_streams_map.end())
      {
        cuda_streams_map[stream_id] = std::move(CreateCudaStream());
      }
      return *cuda_streams_map[stream_id];
    } 
  }  // namespace cuda
}  // namespace micm
