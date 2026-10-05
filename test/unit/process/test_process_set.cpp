#include "test_process_set_policy.hpp"

#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>

using SparseMatrixTest = micm::SparseMatrix<micm::Real>;

using Group1VectorMatrix = micm::VectorMatrix<micm::Real, 1>;
using Group2VectorMatrix = micm::VectorMatrix<micm::Real, 2>;
using Group3VectorMatrix = micm::VectorMatrix<micm::Real, 3>;
using Group4VectorMatrix = micm::VectorMatrix<micm::Real, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<4>>;

TEST(ProcessSet, Matrix)
{
  TestProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest>>();
}

TEST(ProcessSet, VectorMatrix)
{
  TestProcessSet<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestProcessSet<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::ProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestProcessSet<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::ProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestProcessSet<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::ProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(RandomProcessSet, Matrix)
{
  TestRandomSystem<micm::Matrix<micm::Real>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest>>(
      200, 50, 40);
  TestRandomSystem<micm::Matrix<micm::Real>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest>>(
      300, 30, 20);
  TestRandomSystem<micm::Matrix<micm::Real>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest>>(
      400, 100, 80);
}

TEST(ProcessSetAlgebraicVariables, CudaMatrix)
{
  TestAlgebraicMasking<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::ProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::ProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::ProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(ProcessSetFiniteDifferenceJacobian, Matrix)
{
  TestProcessSetFiniteDifferenceJacobian<
      micm::Matrix<micm::Real>,
      SparseMatrixTest,
      micm::ProcessSet<micm::Matrix<micm::Real>, SparseMatrixTest>>();
}