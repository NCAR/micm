// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data
  ///
  /// Data is stored with blocks in the block diagonal matrix as the highest
  /// level structure, then by row, then by non-zero columns in each row.
  class SparseMatrixStandardOrdering
  {
   protected:
    static std::size_t VectorSize(
        std::size_t number_of_blocks,
        const std::vector<std::size_t>& row_ids,
        const std::vector<std::size_t>& row_start)
    {
      return number_of_blocks * row_ids.size();
    };

    static std::size_t VectorIndex(
        std::size_t number_of_blocks,
        const std::vector<std::size_t>& row_ids,
        const std::vector<std::size_t>& row_start,
        std::size_t block,
        std::size_t row,
        std::size_t column)
    {
      if (row >= row_start.size() - 1 || column >= row_start.size() - 1 || block >= number_of_blocks)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      auto begin = std::next(row_ids.begin(), row_start[row]);
      auto end = std::next(row_ids.begin(), row_start[row + 1]);
      auto elem = std::find(begin, end, column);
      if (elem == end)
        throw std::system_error(make_error_code(MicmMatrixErrc::ZeroElementAccess));
      return std::size_t{ (elem - row_ids.begin()) + block * row_ids.size() };
    };

    static void AddToDiagonal(
        const std::vector<std::size_t>& diagonal_ids,
        const std::size_t number_of_blocks,
        const std::size_t block_size,
        auto& data,
        auto value)
    {
      for (std::size_t block_start = 0;
           block_start < number_of_blocks * block_size;
           block_start += block_size)
        for (const auto& i : diagonal_ids)
          data[block_start + i] += value;
    };
  };
}  // namespace micm
