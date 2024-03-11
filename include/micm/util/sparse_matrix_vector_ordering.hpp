// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>

#include "vector_matrix.hpp"

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data to encourage vectorization
  ///
  /// Data is stored with sets of blocks in the block diagonal matrix as the highest
  /// level structure, then by row, then by non-zero columns in each row, then by
  /// individual blocks in the set of blocks.
  ///
  /// The template argument is the number of blocks per set of blocks and should be
  /// approximately the size of the vector register.
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class SparseMatrixVectorOrdering
  {
   protected:
    static std::size_t VectorSize(
        std::size_t number_of_blocks,
        const std::vector<std::size_t>& row_ids,
        const std::vector<std::size_t>& row_start)
    {
      return std::ceil((double)number_of_blocks / (double)L) * L * row_ids.size();
    };

    std::size_t VectorIndex(
        std::size_t number_of_blocks,
        const std::vector<std::size_t>& row_ids,
        const std::vector<std::size_t>& row_start,
        std::size_t block,
        std::size_t row,
        std::size_t column) const
    {
      if (row >= row_start.size() - 1 || column >= row_start.size() - 1 || block >= number_of_blocks)
        throw std::invalid_argument("SparseMatrix element out of range");
      auto begin = std::next(row_ids.begin(), row_start[row]);
      auto end = std::next(row_ids.begin(), row_start[row + 1]);
      auto elem = std::find(begin, end, column);
      if (elem == end)
        throw std::invalid_argument("SparseMatrix zero element access not allowed");
      return std::size_t{ (elem - row_ids.begin()) * L + block % L + (block / L) * L * row_ids.size() };
    };

   public:
    std::size_t GroupVectorSize() const
    {
      return L;
    }

    std::size_t GroupSize(std::size_t number_of_non_zero_elements) const
    {
      return L * number_of_non_zero_elements;
    }

    std::size_t NumberOfGroups(std::size_t number_of_blocks) const
    {
      return std::ceil((double)number_of_blocks / (double)L);
    }
  };

  // Default vectorized SparseMatrix
  template<class T>
  using VectorSparseMatrix = SparseMatrix<T, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

}  // namespace micm