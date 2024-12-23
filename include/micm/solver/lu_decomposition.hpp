// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "lu_decomposition_doolittle.hpp"
#include "lu_decomposition_mozart.hpp"

namespace micm
{
  /// @brief Alias for the default LU decomposition algorithm
  using LuDecomposition = LuDecompositionDoolittle;
}  // namespace micm