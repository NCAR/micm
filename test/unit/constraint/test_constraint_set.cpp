// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "test_constraint_set_policy.hpp"

#include <micm/constraint/constraint_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using namespace micm;
using StandardSparseMatrix = SparseMatrix<micm::Real, SparseMatrixStandardOrdering>;

using Group1VectorMatrix = VectorMatrix<micm::Real, 1>;
using Group2VectorMatrix = VectorMatrix<micm::Real, 2>;
using Group3VectorMatrix = VectorMatrix<micm::Real, 3>;
using Group4VectorMatrix = VectorMatrix<micm::Real, 4>;

using Group1SparseVectorMatrix = SparseMatrix<micm::Real, SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = SparseMatrix<micm::Real, SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = SparseMatrix<micm::Real, SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = SparseMatrix<micm::Real, SparseMatrixVectorOrdering<4>>;

TEST(ConstraintSet, Construction)
{
  TestConstruction<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, ReplaceStateRowsMapsToAlgebraicSpecies)
{
  TestReplaceStateRowsMapsToAlgebraicSpecies<
      Matrix<micm::Real>,
      StandardSparseMatrix,
      ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, NonZeroJacobianElements)
{
  TestNonZeroJacobianElements<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, MultipleConstraints)
{
  TestMultipleConstraints<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, AddForcingTerms)
{
  TestAddForcingTerms<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, SubtractJacobianTerms)
{
  TestSubtractJacobianTerms<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, EmptyConstraintSet)
{
  TestEmptyConstraintSet<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, UnknownSpeciesThrows)
{
  TestUnknownSpeciesThrows<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, ThreeDStateOneConstraint)
{
  TestThreeDStateOneConstraint<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, FourDStateTwoConstraints)
{
  TestFourDStateTwoConstraints<Matrix<micm::Real>, StandardSparseMatrix, ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, CoupledConstraintsSharedSpecies)
{
  TestCoupledConstraintsSharedSpecies<
      Matrix<micm::Real>,
      StandardSparseMatrix,
      ConstraintSet<Matrix<micm::Real>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, VectorizedMatricesRespectGridCellIndexing)
{
  TestVectorizedMatricesRespectGridCellIndexing<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      ConstraintSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix1)
{
  TestConstruction<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestNonZeroJacobianElements<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestAddForcingTerms<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestSubtractJacobianTerms<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix2)
{
  TestConstruction<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestNonZeroJacobianElements<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestAddForcingTerms<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestSubtractJacobianTerms<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix3)
{
  TestConstruction<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestNonZeroJacobianElements<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestAddForcingTerms<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestSubtractJacobianTerms<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}