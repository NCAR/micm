// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/util/cuda_util.cuh>
#include <micm/util/internal_error.hpp>

#include <cuda_runtime.h>
#include <memory>
#include <mutex>
#include <map>
#include <cublas_v2.h>

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

    cublasHandle_t& GetCublasHandle()
    {
      static std::map<int, cublasHandle_t> cublas_handles_map;
      static std::mutex mutex;
      int device_id;
      CHECK_CUDA_ERROR(cudaGetDevice(&device_id), "Failed to get device ID...");
      std::lock_guard<std::mutex> lock(mutex); // lock the mutex and generate a new cublas handle below
      if (auto search = cublas_handles_map.find(device_id); search == cublas_handles_map.end())
      {
        cublasHandle_t handle;
        CHECK_CUBLAS_ERROR(cublasCreate(&handle), "CUBLAS initialization failed..."); // create the cublas handle
        cublas_handles_map[device_id] = handle;                                       // save the cublas handle
      }
      return cublas_handles_map[device_id];
    }
  }  // namespace cuda
}  // namespace micm
