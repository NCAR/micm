// Copyright (C) 2024-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/Process.hpp>
#include <micm/Solver.hpp>
#include <micm/System.hpp>
#include <micm/Util.hpp>
#include <micm/version.hpp>

namespace micm
{
  using DenseMatrixVector = VectorMatrix<double, MICM_DEFAULT_VECTOR_SIZE>;
  using SparseMatrixVector = SparseMatrix<double, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  using DenseMatrixStandard = Matrix<double>;
  using SparseMatrixStandard = SparseMatrix<double, SparseMatrixStandardOrdering>;

  using VectorState = State<DenseMatrixVector, SparseMatrixVector>;
  using StandardState = State<DenseMatrixStandard, SparseMatrixStandard>;

  using RosenbrockVectorType = typename RosenbrockSolverParameters::template SolverType<ProcessSet, LinearSolver<SparseMatrixVector, LuDecomposition>>;
  using Rosenbrock = Solver<RosenbrockVectorType, State<DenseMatrixVector, SparseMatrixVector>>;

  using RosenbrockStandardType = typename RosenbrockSolverParameters::template SolverType<ProcessSet, LinearSolver<SparseMatrixStandard, LuDecomposition>>;
  using RosenbrockStandard = Solver<RosenbrockStandardType, State<DenseMatrixStandard, SparseMatrixStandard>>;

  using BackwardEulerVectorType = typename BackwardEulerSolverParameters::template SolverType<ProcessSet, LinearSolver<SparseMatrixVector, LuDecomposition>>;
  using BackwardEuler = Solver<BackwardEulerVectorType, State<DenseMatrixVector, SparseMatrixVector>>;

  using BackwardEulerStandardType = typename BackwardEulerSolverParameters::template SolverType<ProcessSet, LinearSolver<SparseMatrixStandard, LuDecomposition>>;
  using BackwardEulerStandard = Solver<BackwardEulerStandardType, State<DenseMatrixStandard, SparseMatrixStandard>>;

  using RosenbrockThreeStageBuilder = CpuSolverBuilder<RosenbrockSolverParameters, DenseMatrixVector, SparseMatrixVector>;
      
  using BackwardEulerBuilder = CpuSolverBuilder< BackwardEulerSolverParameters, DenseMatrixVector, SparseMatrixVector, LuDecompositionDoolittle>;
}  // namespace micm