// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <cstddef>
#include <set>
#include <utility>

namespace micm
{
  // annonymous namespace to hide jacobian builder
  template<class SparseMatrixPolicy>
  SparseMatrixPolicy BuildJacobian(
      const std::set<std::pair<Index, Index>>& nonzero_jacobian_elements,
      Index number_of_grid_cells,
      Index state_size,
      bool indexing_only)
  {
    auto builder = SparseMatrixPolicy::Create(state_size).SetNumberOfBlocks(number_of_grid_cells);
    for (const auto& elem : nonzero_jacobian_elements)
    {
      builder = builder.WithElement(elem.first, elem.second);
    }
    // Always include diagonal elements
    for (Index i = 0; i < state_size; ++i)
    {
      builder = builder.WithElement(i, i);
    }
    return SparseMatrixPolicy(builder, indexing_only);
  }
}  // namespace micm
