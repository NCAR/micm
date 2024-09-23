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
  testDenseMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<SparseMatrixTest>>();
}

TEST(LinearSolver, RandomMatrixStandardOrdering)
{
  testRandomMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<SparseMatrixTest>>(5);
}

TEST(LinearSolver, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<SparseMatrixTest>>(5);
}

TEST(LinearSolver, DiagonalMarkowitzReorder)
{
  testMarkowitzReordering<micm::Matrix<int>, SparseMatrixTest>();
}

TEST(LinearSolver, StandardOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for(auto initial_value : initial_values)
    testExtremeInitialValue<DenseMatrixTest, SparseMatrixTest, micm::LinearSolver<SparseMatrixTest>>(5, initial_value);
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
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<Group1SparseVectorMatrix>>();
  testDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<Group2SparseVectorMatrix>>();
  testDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<Group3SparseVectorMatrix>>();
  testDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<Group4SparseVectorMatrix>>();
}

TEST(LinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<Group1SparseVectorMatrix>>(5);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<Group2SparseVectorMatrix>>(5);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<Group3SparseVectorMatrix>>(5);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, VectorOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for(auto initial_value : initial_values) {
    testExtremeInitialValue<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<Group1SparseVectorMatrix>>(1, initial_value);
    testExtremeInitialValue<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<Group2SparseVectorMatrix>>(2, initial_value);
    testExtremeInitialValue<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<Group3SparseVectorMatrix>>(5, initial_value);
    testExtremeInitialValue<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<Group4SparseVectorMatrix>>(5, initial_value);
  }
}

TEST(LinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<Group1SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<Group2SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<Group3SparseVectorMatrix>>(5);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolver, VectorDiagonalMarkowitzReordering)
{
  testMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
