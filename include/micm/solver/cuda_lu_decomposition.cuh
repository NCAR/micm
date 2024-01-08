// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform LU decomposition on the device
    std::chrono::nanoseconds DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            LuDecomposeConst* devptr);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecomposition" to the device;
    /// Note that if we want to allocate "devptr" inside this function,
    ///   it must be declared as "LuDecomposeConstDevice*&"; otherwise,
    ///   passing "devptr->d_niLU_" as an argument to the CUDA kernel
    ///   will trigger a segmentation fault;
    void CopyConstData(LuDecomposeConst* hostptr, LuDecomposeConst*& devptr);
   
    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeConst*& devptr);

  }  // end of namespace cuda
}    // end of namespace micm