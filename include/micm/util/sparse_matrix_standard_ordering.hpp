// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "sparse_matrix_standard_ordering_compressed_sparse_column.hpp"
#include "sparse_matrix_standard_ordering_compressed_sparse_row.hpp"
#include "view_category.hpp"

namespace micm
{

  /// @brief Alias for the default sparse matrix standard ordering
  using SparseMatrixStandardOrdering = SparseMatrixStandardOrderingCompressedSparseRow;

  template<class T, class OrderingPolicy>
  class SparseMatrix;

  /// @brief Standard ordering row sparse matrices always use simple grouping (L==1)
  template<typename T>
  struct GroupingStrategy<SparseMatrix<SparseMatrixStandardOrderingCompressedSparseRow, T>>
  {
    using type = SimpleGroupingTag;
  };

  /// @brief Standard ordering column sparse matrices always use simple grouping (L==1)
  template<typename T>
  struct GroupingStrategy<SparseMatrix<SparseMatrixStandardOrderingCompressedSparseColumn, T>>
  {
    using type = SimpleGroupingTag;
  };

}  // namespace micm
