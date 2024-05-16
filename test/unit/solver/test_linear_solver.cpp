#include "test_linear_solver_policy.hpp"

#include <micm/solver/linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <functional>

using FloatingPointType = double;

using DenseMatrixTest = micm::Matrix<FloatingPointType>;
using SparseMatrixTest = micm::SparseMatrix<FloatingPointType>;

TEST(LinearSolver, DenseMatrixStandardOrdering)
{
  testDenseMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<FloatingPointType, SparseMatrixTest>>();
}

TEST(LinearSolver, RandomMatrixStandardOrdering)
{
  testRandomMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<FloatingPointType, SparseMatrixTest>>(5);
}

TEST(LinearSolver, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<FloatingPointType, SparseMatrixTest>>(5);
}

TEST(LinearSolver, DiagonalMarkowitzReorder)
{
  testMarkowitzReordering<micm::Matrix<int>, SparseMatrixTest>();
}

using Group1VectorMatrix = micm::VectorMatrix<FloatingPointType, 1>;
using Group2VectorMatrix = micm::VectorMatrix<FloatingPointType, 2>;
using Group3VectorMatrix = micm::VectorMatrix<FloatingPointType, 3>;
using Group4VectorMatrix = micm::VectorMatrix<FloatingPointType, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<4>>;

TEST(LinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group1SparseVectorMatrix>>();
  testDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group2SparseVectorMatrix>>();
  testDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group3SparseVectorMatrix>>();
  testDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group4SparseVectorMatrix>>();
}

TEST(LinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group1SparseVectorMatrix>>(5);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group2SparseVectorMatrix>>(5);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group3SparseVectorMatrix>>(5);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group1SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group2SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group3SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<FloatingPointType, Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, VectorDiagonalMarkowitzReordering)
{
  testMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
