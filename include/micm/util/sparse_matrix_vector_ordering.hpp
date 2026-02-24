// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "sparse_matrix_vector_ordering_compressed_sparse_column.hpp"
#include "sparse_matrix_vector_ordering_compressed_sparse_row.hpp"
#include "view_category.hpp"

namespace micm
{
  /// @brief Alias for the default sparse matrix vector ordering
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using SparseMatrixVectorOrdering = SparseMatrixVectorOrderingCompressedSparseRow<L>;

  // Default vectorized SparseMatrix
  using DefaultVectorSparseMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  // Define the following two functions that only work for the CudaMatrix; the if constexpr statement is evaluated at
  // compile-time Reference: https://www.modernescpp.com/index.php/using-requires-expression-in-c-20-as-a-standalone-feature/
  template<class MatrixPolicy>
  void CheckCopyToDevice(MatrixPolicy& matrix)
  {
    if constexpr (requires {
                    { matrix.CopyToDevice() } -> std::same_as<void>;
                  })
      matrix.CopyToDevice();
  }

  template<class MatrixPolicy>
  void CheckCopyToHost(MatrixPolicy& matrix)
  {
    if constexpr (requires {
                    { matrix.CopyToHost() } -> std::same_as<void>;
                  })
      matrix.CopyToHost();
  }

  // ============================================================================
  // Grouping Strategy Specializations for Vector Ordering
  // ============================================================================

  /// @brief Vector ordering row sparse matrices use simple grouping when L==1, tiered when L>1
  template<std::size_t L, typename T>
  struct GroupingStrategy<SparseMatrix<SparseMatrixVectorOrderingCompressedSparseRow<L>, T>>
  {
    using type = std::conditional_t<L == 1, SimpleGroupingTag, TieredGroupingTag>;
  };

  /// @brief Vector ordering column sparse matrices use simple grouping when L==1, tiered when L>1
  template<std::size_t L, typename T>
  struct GroupingStrategy<SparseMatrix<SparseMatrixVectorOrderingCompressedSparseColumn<L>, T>>
  {
    using type = std::conditional_t<L == 1, SimpleGroupingTag, TieredGroupingTag>;
  };

}  // namespace micm