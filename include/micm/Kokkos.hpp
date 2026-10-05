// Copyright (C) 2024-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/CPU.hpp>
#include <micm/kokkos/solver/kokkos_solver_builder.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/solver/state.hpp>

namespace micm
{
  using KokkosDenseReal = KokkosDenseMatrix<Real, MICM_DEFAULT_VECTOR_SIZE>;
  using KokkosSparseReal = KokkosSparseMatrix<Real, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  using KokkosState = State<KokkosDenseReal, KokkosSparseReal, LuDecompositionMozartInPlace<KokkosSparseReal>>;

  using KokkosRosenbrockType = typename RosenbrockSolverParameters::template SolverType<
      ProcessSet<KokkosDenseReal, KokkosSparseReal>,
      LinearSolverInPlace<KokkosDenseReal, KokkosSparseReal, LuDecompositionMozartInPlace<KokkosSparseReal>>,
      ConstraintSet<KokkosDenseReal, KokkosSparseReal>>;
  using KokkosRosenbrock = Solver<KokkosRosenbrockType, KokkosState>;

  using KokkosBackwardEulerType = typename BackwardEulerSolverParameters::template SolverType<
      ProcessSet<KokkosDenseReal, KokkosSparseReal>,
      LinearSolverInPlace<KokkosDenseReal, KokkosSparseReal, LuDecompositionMozartInPlace<KokkosSparseReal>>,
      ConstraintSet<KokkosDenseReal, KokkosSparseReal>>;
  using KokkosBackwardEuler = Solver<KokkosBackwardEulerType, KokkosState>;

}  // namespace micm
