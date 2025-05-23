// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "jit_lu_decomposition_doolittle.hpp"

namespace micm
{
  /// @brief Alias for the default JIT LU decomposition algorithm
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using JitLuDecomposition = JitLuDecompositionDoolittle<L>;
}  // namespace micm