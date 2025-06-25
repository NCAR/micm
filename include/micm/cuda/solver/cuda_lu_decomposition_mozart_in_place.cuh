// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform LU decomposition on the device
    void DecomposeKernelDriver(CudaMatrixParam& ALU_param, const LuDecomposeMozartInPlaceParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLuDecompositionMozartInPlace" to the device;
    LuDecomposeMozartInPlaceParam CopyConstData(LuDecomposeMozartInPlaceParam& hoststruct);

    /// This is the function that will delete the constant data
    ///   members of class "CudaLuDecompositionMozartInPlace" on the device
    void FreeConstData(LuDecomposeMozartInPlaceParam& devstruct);

  }  // end of namespace cuda
}  // end of namespace micm
