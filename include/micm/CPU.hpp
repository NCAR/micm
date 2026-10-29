// Copyright (C) 2024-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/Process.hpp>
#include <micm/Solver.hpp>
#include <micm/System.hpp>
#include <micm/Util.hpp>
#include <micm/external_model.hpp>
#include <micm/util/types.hpp>
#include <micm/version.hpp>

namespace micm
{
  using DenseMatrixVector = VectorMatrix<Real, MICM_DEFAULT_VECTOR_SIZE>;
  using SparseMatrixVector = SparseMatrix<Real, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  using DenseMatrixStandard = Matrix<Real>;
  using SparseMatrixStandard = SparseMatrix<Real, SparseMatrixStandardOrdering>;

  using LuDecompositionVector = LuDecompositionMozartInPlace<SparseMatrixVector>;
  using LinearSolverVector = LinearSolverInPlace<DenseMatrixVector, SparseMatrixVector, LuDecompositionVector>;

  using VectorState = State<DenseMatrixVector, SparseMatrixVector, LuDecompositionVector>;
  using LuDecompositionStandard = LuDecompositionInPlace<SparseMatrixStandard>;
  using LinearSolverStandard = LinearSolverInPlace<DenseMatrixStandard, SparseMatrixStandard, LuDecompositionStandard>;

  using StandardState = State<DenseMatrixStandard, SparseMatrixStandard, LuDecompositionStandard>;

  using RosenbrockVectorType = typename RosenbrockSolverParameters::template SolverType<
      ProcessSet<DenseMatrixVector, SparseMatrixVector>,
      LinearSolverVector,
      ConstraintSet<DenseMatrixVector, SparseMatrixVector>>;
  using Rosenbrock = Solver<RosenbrockVectorType, VectorState>;

  using RosenbrockStandardType = typename RosenbrockSolverParameters::template SolverType<
      ProcessSet<DenseMatrixStandard, SparseMatrixStandard>,
      LinearSolverStandard,
      ConstraintSet<DenseMatrixStandard, SparseMatrixStandard>>;
  using RosenbrockStandard = Solver<RosenbrockStandardType, StandardState>;

  using BackwardEulerVectorType = typename BackwardEulerSolverParameters::template SolverType<
      ProcessSet<DenseMatrixVector, SparseMatrixVector>,
      LinearSolverVector,
      ConstraintSet<DenseMatrixVector, SparseMatrixVector>>;
  using BackwardEuler = Solver<BackwardEulerVectorType, VectorState>;

  using BackwardEulerStandardType = typename BackwardEulerSolverParameters::template SolverType<
      ProcessSet<DenseMatrixStandard, SparseMatrixStandard>,
      LinearSolverStandard,
      ConstraintSet<DenseMatrixStandard, SparseMatrixStandard>>;
  using BackwardEulerStandard = Solver<BackwardEulerStandardType, StandardState>;

  using RosenbrockThreeStageBuilder =
      CpuSolverBuilder<RosenbrockSolverParameters, DenseMatrixVector, SparseMatrixVector, LuDecompositionVector>;
  using BackwardEulerBuilder =
      CpuSolverBuilder<BackwardEulerSolverParameters, DenseMatrixVector, SparseMatrixVector, LuDecompositionVector>;
}  // namespace micm