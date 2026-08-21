// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/solver/linear_solver_in_place.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

namespace micm
{
  /// @brief Builder of Kokkos-backed solvers
  /// @tparam SolverParametersPolicy Policy for the solver parameters struct
  /// @tparam L Vector dimension
  template<class SolverParametersPolicy, Index L = MICM_KOKKOS_DEFAULT_VECTOR_SIZE>
  using KokkosSolverBuilder = SolverBuilder<
      SolverParametersPolicy,
      KokkosDenseMatrix<Real, L>,
      KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>,
      ProcessSet<KokkosDenseMatrix<Real, L>, KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>>,
      LuDecompositionMozartInPlace<KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>>,
      LinearSolverInPlace<
          KokkosDenseMatrix<Real, L>,
          KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>,
          LuDecompositionMozartInPlace<KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>>>,
      State<
          KokkosDenseMatrix<Real, L>,
          KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>,
          LuDecompositionMozartInPlace<KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<L>>>>>;
}  // namespace micm