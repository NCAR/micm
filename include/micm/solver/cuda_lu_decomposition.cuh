// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm
{
  namespace cuda
  {
    void DecomposeKernelDriver(
            CudaSparseMatrixParam& sparseMatrix, 
            CudaSolverParam& solver);
  }  // namespace cuda
}  // namespace micm
