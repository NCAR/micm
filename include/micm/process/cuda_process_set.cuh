// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    std::chrono::nanoseconds AddForcingTermsKernelDriver(
        CudaMatrixParam& matrixParam,
        CudaProcessSetParam& processSet); 

    std::chrono::nanoseconds AddJacobianTermsKernelDriver(
        CudaMatrixParam& matrixParam, 
        CudaSparseMatrixParam& sparseMatrix, 
        CudaProcessSetParam& processSet);
  }  // namespace cuda
}  // namespace micm
