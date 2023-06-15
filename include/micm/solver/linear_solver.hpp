// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief A general-use sparse-matrix linear solver
  class LinearSolver
  {
    public:
    /// @brief default constructor
    LinearSolver();

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    LinearSolver(const SparseMatrix& matrix);


  };
} // namespace micm