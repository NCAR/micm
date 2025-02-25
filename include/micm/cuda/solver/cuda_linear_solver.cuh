// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the host function that will call the CUDA kernel
    ///   to perform the "solve" function on the device
    void SolveKernelDriver(
        CudaMatrixParam& x_param,
        const CudaMatrixParam& L_param,
        const CudaMatrixParam& U_param,
        const LinearSolverParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLinearSolver" to the device;
    LinearSolverParam CopyConstData(LinearSolverParam& hoststruct);

    /// This is the function that will delete the constant data
    ///   members of class "CudaLinearSolver" on the device
    void FreeConstData(LinearSolverParam& devstruct);
  }  // namespace cuda
}  // namespace micm
