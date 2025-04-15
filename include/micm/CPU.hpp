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
  using DenseMatrixVector = micm::VectorMatrix<double, MICM_DEFAULT_VECTOR_SIZE>;
  using SparseMatrixVector = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  using DenseMatrixStandard = micm::Matrix<double>;
  using SparseMatrixStandard = micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>;

  using VectorState = micm::State<DenseMatrixVector, SparseMatrixVector>;
  using StandardState = micm::State<DenseMatrixStandard, SparseMatrixStandard>;

  using RosenbrockVectorType = typename micm::RosenbrockSolverParameters::
      template SolverType<micm::ProcessSet, micm::LinearSolver<SparseMatrixVector, micm::LuDecomposition>>;
  using Rosenbrock = micm::Solver<RosenbrockVectorType, micm::State<DenseMatrixVector, SparseMatrixVector>>;

  using RosenbrockStandardType = typename micm::RosenbrockSolverParameters::
      template SolverType<micm::ProcessSet, micm::LinearSolver<SparseMatrixStandard, micm::LuDecomposition>>;
  using RosenbrockStandard = micm::Solver<RosenbrockStandardType, micm::State<DenseMatrixStandard, SparseMatrixStandard>>;

  using BackwardEulerVectorType = typename micm::BackwardEulerSolverParameters::
      template SolverType<micm::ProcessSet, micm::LinearSolver<SparseMatrixVector, micm::LuDecomposition>>;
  using BackwardEuler = micm::Solver<BackwardEulerVectorType, micm::State<DenseMatrixVector, SparseMatrixVector>>;

  using BackwardEulerStandardType = typename micm::BackwardEulerSolverParameters::
      template SolverType<micm::ProcessSet, micm::LinearSolver<SparseMatrixStandard, micm::LuDecomposition>>;
  using BackwardEulerStandard =
      micm::Solver<BackwardEulerStandardType, micm::State<DenseMatrixStandard, SparseMatrixStandard>>;

  using RosenbrockThreeStageBuilder =
      micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DenseMatrixVector, SparseMatrixVector>;
  using BackwardEulerBuilder = micm::CpuSolverBuilder<
      micm::BackwardEulerSolverParameters,
      DenseMatrixVector,
      SparseMatrixVector,
      micm::LuDecompositionDoolittle>;
}  // namespace micm