// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "sparse_matrix_standard_ordering_compressed_sparse_column.hpp"
#include "sparse_matrix_standard_ordering_compressed_sparse_row.hpp"

namespace micm
{

  /// @brief Alias for the default sparse matrix standard ordering
  using SparseMatrixStandardOrdering = SparseMatrixStandardOrderingCompressedSparseRow;

}  // namespace micm
