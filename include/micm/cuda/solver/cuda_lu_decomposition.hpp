// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "cuda_lu_decomposition_mozart_in_place.hpp"

namespace micm
{
  /// @brief Alias for the default CUDA LU decomposition algorithm
  using CudaLuDecomposition = CudaLuDecompositionMozartInPlace;
}  // namespace micm
