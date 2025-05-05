// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "sparse_matrix_vector_ordering_compressed_sparse_column.hpp"
#include "sparse_matrix_vector_ordering_compressed_sparse_row.hpp"

namespace micm
{

  /// @brief Alias for the default sparse matrix vector ordering
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using SparseMatrixVectorOrdering = SparseMatrixVectorOrderingCompressedSparseRow<L>;

  // Default vectorized SparseMatrix
  using DefaultVectorSparseMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;
}  // namespace micm