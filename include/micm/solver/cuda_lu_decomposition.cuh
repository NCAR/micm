// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>
namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform LU decomposition on the device
    std::chrono::nanoseconds DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            const LuDecomposeParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecomposition" to the device;
    LuDecomposeParam CopyConstData(LuDecomposeParam& hoststruct);
   
    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeParam& devstruct);

  }  // end of namespace cuda
}    // end of namespace micm
