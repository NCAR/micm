// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace micm
{
  // annonymous namespace to hide jacobian builder
  template<template<class> class SparseMatrixPolicy>
  SparseMatrixPolicy<double> build_jacobian(
      std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements,
      size_t number_of_grid_cells,
      size_t state_size)
  {
    auto builder = SparseMatrixPolicy<double>::create(state_size).number_of_blocks(number_of_grid_cells);
    for (auto& elem : nonzero_jacobian_elements)
      builder = builder.with_element(elem.first, elem.second);
    // Always include diagonal elements
    for (std::size_t i = 0; i < state_size; ++i)
      builder = builder.with_element(i, i);

    return SparseMatrixPolicy<double>(builder);
  }
}  // namespace micm