// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "solver_builder.hpp"
#include <micm/solver/jit_linear_solver.hpp>
#include <micm/solver/jit_lu_decomposition.hpp>
#include <micm/process/jit_process_set.hpp>

namespace micm
{
  /// @brief Builder of CPU-based solvers optimized for a particular chemical system using JIT compilation
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam T Primary data type
  /// @tparam L Vector size
  ///
  /// JIT-compiled solvers only work with vector-ordered matrices
  template<class SolverParametersPolicy, class T = double, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using JitSolverBuilder = SolverBuilder<SolverParametersPolicy, VectorMatrix<T, L>, SparseMatrix<T, SparseMatrixVectorOrdering<L>>, JitProcessSet, JitLinearSolver<SparseMatrixVectorOrdering<L>>, JitLuDecomposition>>;
} // namespace micm