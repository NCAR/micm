/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
 #include <micm/util/cuda_util.cuh>
 #include <micm/util/internal_error.hpp>
 #include <cuda_runtime.h>
 
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
   }  // namespace cuda
 }  // namespace micm
 