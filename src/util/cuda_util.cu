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
  }  // namespace cuda
}  // namespace micm
