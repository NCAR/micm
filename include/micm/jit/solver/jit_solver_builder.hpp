/* Copyright (C) 2023-2025 University Corporation for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/jit/process/jit_process_set.hpp>
#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/jit/solver/jit_lu_decomposition.hpp>
#include <micm/solver/solver_builder.hpp>

namespace micm
{
  /// @brief Builder of CPU-based solvers optimized for a particular chemical system using JIT compilation
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam L Vector size
  ///
  /// JIT-compiled solvers only work with vector-ordered matrices
  template<class SolverParametersPolicy, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using JitSolverBuilder = SolverBuilder<
      SolverParametersPolicy,
      VectorMatrix<double, L>,
      SparseMatrix<double, SparseMatrixVectorOrdering<L>>,
      JitProcessSet<L>,
      JitLuDecomposition<L>,
      JitLinearSolver<L, SparseMatrix<double, SparseMatrixVectorOrdering<L>>, JitLuDecomposition<L>>,
      State<VectorMatrix<double, L>, SparseMatrix<double, SparseMatrixVectorOrdering<L>>, JitLuDecomposition<L>>>;
}  // namespace micm