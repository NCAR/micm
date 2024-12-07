// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "cuda_lu_decomposition_doolittle.hpp"

namespace micm
{
  /// @brief Alias for the default CUDA LU decomposition algorithm
  using CudaLuDecomposition = CudaLuDecompositionDoolittle;
}  // namespace micm
