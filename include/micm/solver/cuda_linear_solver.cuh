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
    ///   to perform the "solve" function on the device
    std::chrono::nanoseconds SolveKernelDriver(
           CudaSparseMatrixParam& sparseMatrix, 
           CudaMatrixParam& denseMatrix,
           const LinearSolverParam& devstruct);

    /// This is the function that will copy the constant data
    ///   members of class "CudaLinearSolver" to the device;
    LinearSolverParam CopyConstData(LinearSolverParam& hoststruct);
   
    /// This is the function that will delete the constant data
    ///   members of class "CudaLinearSolver" on the device
    void FreeConstData(LinearSolverParam& devstruct);
  }  // namespace cuda
}  // namespace micm
