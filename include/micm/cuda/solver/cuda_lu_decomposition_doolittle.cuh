// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform LU decomposition on the device
    void DecomposeKernelDriver(
        const CudaMatrixParam& A_param,
        CudaMatrixParam& L_param,
        CudaMatrixParam& U_param,
        const LuDecomposeDoolittleParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecomposition" to the device;
    LuDecomposeDoolittleParam CopyConstData(LuDecomposeDoolittleParam& hoststruct);

    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecomposition" on the device
    void FreeConstData(LuDecomposeDoolittleParam& devstruct);

  }  // end of namespace cuda
}  // end of namespace micm
