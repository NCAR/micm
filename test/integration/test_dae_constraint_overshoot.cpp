// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "dae_constraint_overshoot_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using DenseMatrix = micm::VectorMatrix<micm::Real, MICM_DEFAULT_VECTOR_SIZE>;
using StdSparseMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;
using Builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DenseMatrix, StdSparseMatrix>;

TEST(DAEConstraintOvershoot, AlgebraicVariableStaysNonNegative)
{
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  TestAlgebraicVariableStaysNonNegative(Builder(options));
}

TEST(DAEConstraintOvershoot, EquilibriumPlusConservation)
{
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  TestEquilibriumPlusConservation(Builder(options));
}

TEST(DAEConstraintOvershoot, AllRosenbrockOrdersConstrained)
{
  TestAllRosenbrockOrdersConstrained([](const micm::RosenbrockSolverParameters& options) { return Builder(options); });
}
