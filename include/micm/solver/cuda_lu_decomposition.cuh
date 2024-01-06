// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  /// This is a forward declaration of class "CudaLuDecomposion";
  /// The linker will know where the acutual class is defined and used;
  class CudaLuDecomposition;
  namespace cuda
  {
    /// This is the CPU function that calls the
    /// CUDA kernel to perform LU decomposition;
    std::chrono::nanoseconds DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            LuDecomposeConstDevice* devptr);

    /// This is the function to copy the constant data members 
    /// of objects with the "CudaLuDecomposition" type to the device
    void CopyConstData(LuDecomposeConstHost* hostptr, LuDecomposeConstDevice* devptr);
   
    /// This is the function to delete the constant data members 
    /// of objects with the "CudaLuDecomposition" type on the device
    void FreeConstData(LuDecomposeConstDevice* devptr);

  }  // end of namespace cuda
}    // end of namespace micm
