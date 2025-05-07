// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform the "solve" function on the device
    void
    SolveKernelDriver(CudaMatrixParam& x_param, const CudaMatrixParam& ALU_param, const LinearSolverInPlaceParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLinearSolverInPlace" to the device;
    LinearSolverInPlaceParam CopyConstData(LinearSolverInPlaceParam& hoststruct);

    /// This is the function that will delete the constant data
    ///   members of class "CudaLinearSolverInPlace" on the device
    void FreeConstData(LinearSolverInPlaceParam& devstruct);
  }  // namespace cuda
}  // namespace micm
