// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../constraint/test_constraint_set_policy.hpp"

#include <micm/constraint/constraint_set.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using namespace micm;
using Group1VectorMatrix = KokkosDenseMatrix<micm::Real, 1>;
using Group2VectorMatrix = KokkosDenseMatrix<micm::Real, 2>;
using Group3VectorMatrix = KokkosDenseMatrix<micm::Real, 3>;
using Group4VectorMatrix = KokkosDenseMatrix<micm::Real, 4>;

using Group1SparseVectorMatrix = KokkosSparseMatrix<micm::Real, SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = KokkosSparseMatrix<micm::Real, SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = KokkosSparseMatrix<micm::Real, SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = KokkosSparseMatrix<micm::Real, SparseMatrixVectorOrdering<4>>;

TEST(ConstraintSet, Construction)
{
  TestConstruction<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, ReplaceStateRowsMapsToAlgebraicSpecies)
{
  TestReplaceStateRowsMapsToAlgebraicSpecies<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, NonZeroJacobianElements)
{
  TestNonZeroJacobianElements<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, MultipleConstraints)
{
  TestMultipleConstraints<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, AddForcingTerms)
{
  TestAddForcingTerms<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, SubtractJacobianTerms)
{
  TestSubtractJacobianTerms<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, EmptyConstraintSet)
{
  TestEmptyConstraintSet<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, UnknownSpeciesThrows)
{
  TestUnknownSpeciesThrows<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, ThreeDStateOneConstraint)
{
  TestThreeDStateOneConstraint<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, FourDStateTwoConstraints)
{
  TestFourDStateTwoConstraints<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
}

TEST(ConstraintSet, CoupledConstraintsSharedSpecies)
{
  TestCoupledConstraintsSharedSpecies<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      ConstraintSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
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

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}