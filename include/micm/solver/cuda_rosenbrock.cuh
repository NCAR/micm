// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <micm/util/cuda_param.hpp>
#include <vector>
namespace micm{
    namespace cuda{
  std::chrono::nanoseconds AlphaMinusJacobianDriver(
                        CudaSparseMatrixParam& sparseMatrix,
                        const std::vector<size_t> jacobian_diagonal_elements,
                        double alpha);

    }//end cuda
}//end micm