// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "linear_constraint_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using namespace micm;
using DenseMatrix = VectorMatrix<Real, MICM_DEFAULT_VECTOR_SIZE>;
using StdSparseMatrix = SparseMatrix<Real, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

TEST(DAESolveWithConstraint, TerminatorAndRobertson)
{
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  TestTerminatorAndRobertson(
      micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DenseMatrix, StdSparseMatrix>(options));
}
