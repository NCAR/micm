#include "test_linear_solver_policy.hpp"

#include <micm/solver/linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <functional>

using DenseMatrixTest = micm::Matrix<double>;
using SparseMatrixTest = micm::SparseMatrix<double>;

TEST(LinearSolver, DenseMatrixStandardOrdering)
{
  testDenseMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<double, micm::SparseMatrix<>>>();
}

TEST(LinearSolver, RandomMatrixStandardOrdering)
{
  testRandomMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<double, micm::SparseMatrix<>>>(5);
}

TEST(LinearSolver, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<double, micm::SparseMatrix<>>>(5);
}

TEST(LinearSolver, DiagonalMarkowitzReorder)
{
  testMarkowitzReordering<micm::Matrix<int>, SparseMatrixTest>();
}

using Group1VectorMatrix = micm::VectorMatrix<double, 1>;
using Group2VectorMatrix = micm::VectorMatrix<double, 2>;
using Group3VectorMatrix = micm::VectorMatrix<double, 3>;
using Group4VectorMatrix = micm::VectorMatrix<double, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(LinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>();
  testDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>();
  testDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>();
  testDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>();
}

TEST(LinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(5);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(5);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(5);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, VectorDiagonalMarkowitzReordering)
{
  testMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
