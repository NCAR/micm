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
    /// Note that referencing "devstruct" as "LuDecomposeConst&"
    ///   will cause an error like `error: qualifiers dropped in 
    ///   binding reference of type "LuDecomposeConst &" to 
    ///   initializer of type "const LuDecomposeConst"`
    std::chrono::nanoseconds DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            LuDecomposeConst devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecomposition" to the device;
    LuDecomposeConst CopyConstData(LuDecomposeConst& hoststruct);
   
    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeConst& devstruct);

  }  // end of namespace cuda
}    // end of namespace micm
