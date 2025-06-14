// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>

#include <cstddef>
#include <set>
#include <utility>

namespace micm
{
  // annonymous namespace to hide jacobian builder
  template<class SparseMatrixPolicy>
  SparseMatrixPolicy BuildJacobian(
      const std::set<std::pair<std::size_t, std::size_t>>& nonzero_jacobian_elements,
      std::size_t number_of_grid_cells,
      std::size_t state_size,
      bool indexing_only)
  {
    MICM_PROFILE_FUNCTION();

    auto builder = SparseMatrixPolicy::Create(state_size).SetNumberOfBlocks(number_of_grid_cells);
    for (auto& elem : nonzero_jacobian_elements)
      builder = builder.WithElement(elem.first, elem.second);
    // Always include diagonal elements
    for (std::size_t i = 0; i < state_size; ++i)
    {
      builder = builder.WithElement(i, i);
    }
    return SparseMatrixPolicy(builder, indexing_only);
  }
}  // namespace micm
