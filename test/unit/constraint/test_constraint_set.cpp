// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "test_constraint_set_policy.hpp"

#include <micm/constraint/constraint_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using namespace micm;
using StandardSparseMatrix = SparseMatrix<double, SparseMatrixStandardOrdering>;

using Group1VectorMatrix = VectorMatrix<double, 1>;
using Group2VectorMatrix = VectorMatrix<double, 2>;
using Group3VectorMatrix = VectorMatrix<double, 3>;
using Group4VectorMatrix = VectorMatrix<double, 4>;

using Group1SparseVectorMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = SparseMatrix<double, SparseMatrixVectorOrdering<4>>;

TEST(ConstraintSet, Construction)
{
  testConstruction<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, ReplaceStateRowsMapsToAlgebraicSpecies)
{
  testReplaceStateRowsMapsToAlgebraicSpecies<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, NonZeroJacobianElements)
{
  testNonZeroJacobianElements<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, MultipleConstraints)
{
  testMultipleConstraints<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, AddForcingTerms)
{
  testAddForcingTerms<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, SubtractJacobianTerms)
{
  testSubtractJacobianTerms<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, EmptyConstraintSet)
{
  testEmptyConstraintSet<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, UnknownSpeciesThrows)
{
  testUnknownSpeciesThrows<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, ThreeDStateOneConstraint)
{
  testThreeDStateOneConstraint<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, FourDStateTwoConstraints)
{
  testFourDStateTwoConstraints<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, CoupledConstraintsSharedSpecies)
{
  testCoupledConstraintsSharedSpecies<Matrix<double>, StandardSparseMatrix, ConstraintSet<Matrix<double>, StandardSparseMatrix>>();
}

TEST(ConstraintSet, VectorizedMatricesRespectGridCellIndexing)
{
  testVectorizedMatricesRespectGridCellIndexing<Group4VectorMatrix, Group4SparseVectorMatrix, ConstraintSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix1)
{
  testConstruction<Group1VectorMatrix, Group1SparseVectorMatrix, ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  testNonZeroJacobianElements<Group1VectorMatrix, Group1SparseVectorMatrix, ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  testAddForcingTerms<Group1VectorMatrix, Group1SparseVectorMatrix, ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  testSubtractJacobianTerms<Group1VectorMatrix, Group1SparseVectorMatrix, ConstraintSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix2)
{
  testConstruction<Group2VectorMatrix, Group2SparseVectorMatrix, ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  testNonZeroJacobianElements<Group2VectorMatrix, Group2SparseVectorMatrix, ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  testAddForcingTerms<Group2VectorMatrix, Group2SparseVectorMatrix, ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  testSubtractJacobianTerms<Group2VectorMatrix, Group2SparseVectorMatrix, ConstraintSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
}

TEST(ConstraintSet, VectorMatrix3)
{
  testConstruction<Group3VectorMatrix, Group3SparseVectorMatrix, ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  testNonZeroJacobianElements<Group3VectorMatrix, Group3SparseVectorMatrix, ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  testAddForcingTerms<Group3VectorMatrix, Group3SparseVectorMatrix, ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  testSubtractJacobianTerms<Group3VectorMatrix, Group3SparseVectorMatrix, ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}