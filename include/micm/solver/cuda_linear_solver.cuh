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
    std::chrono::nanoseconds SolveKernelDriver(
     CudaLinearSolverParam& linearSolver,
     CudaSparseMatrixParam& sparseMatrix, 
     CudaMatrixParam& denseMatrix);
  }  // namespace cuda
}  // namespace micm
