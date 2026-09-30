// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "equilibrium_constraint_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

namespace
{
  using DenseMatrix = micm::VectorMatrix<micm::Real, MICM_DEFAULT_VECTOR_SIZE>;
  using StdSparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  auto MakeBuilder(const micm::RosenbrockSolverParameters& options)
  {
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DenseMatrix, StdSparseMatrix>(options);
  }
}  // namespace

TEST(EquilibriumIntegration, SetConstraintsAPIWorks)
{
  TestSetConstraintsAPIWorks(MakeBuilder);
}

TEST(EquilibriumIntegration, SetConstraintsAPIMultipleConstraints)
{
  TestSetConstraintsAPIMultipleConstraints(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithConstraint)
{
  TestDAESolveWithConstraint(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithConstraintAndReorderState)
{
  TestDAESolveWithConstraintAndReorderState(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithTwoCoupledConstraints)
{
  TestDAESolveWithTwoCoupledConstraints(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithNonUnitStoichiometry)
{
  TestDAESolveWithNonUnitStoichiometry(MakeBuilder);
}
