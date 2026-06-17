#include "test_linear_solver_in_place_policy.hpp"

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

TEST(LinearSolverInPlace, DenseMatrixStandardOrdering)
{
  TestDenseMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolverInPlace<SparseMatrixTest>>();
}

TEST(LinearSolverInPlace, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolverInPlace<SparseMatrixTest>>(5);
}

TEST(LinearSolverInPlace, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<DenseMatrixTest, SparseMatrixTest, micm::LinearSolverInPlace<SparseMatrixTest>>(5);
}

TEST(LinearSolverInPlace, DiagonalMarkowitzReorder)
{
  TestMarkowitzReordering<micm::Matrix<int>, SparseMatrixTest>();
}

TEST(LinearSolverInPlace, StandardOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto initial_value : initial_values)
  {
    TestExtremeInitialValue<DenseMatrixTest, SparseMatrixTest, micm::LinearSolverInPlace<SparseMatrixTest>>(
        5, initial_value);
  }
}

using Group1VectorMatrix = micm::VectorMatrix<FloatingPointType, 1>;
using Group2VectorMatrix = micm::VectorMatrix<FloatingPointType, 2>;
using Group3VectorMatrix = micm::VectorMatrix<FloatingPointType, 3>;
using Group4VectorMatrix = micm::VectorMatrix<FloatingPointType, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<4>>;

TEST(LinearSolverInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolverInPlace<Group1SparseVectorMatrix>>();
  TestDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolverInPlace<Group2SparseVectorMatrix>>();
  TestDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolverInPlace<Group3SparseVectorMatrix>>();
  TestDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolverInPlace<Group4SparseVectorMatrix>>();
}

TEST(LinearSolverInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolverInPlace<Group1SparseVectorMatrix>>(5);
  TestRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolverInPlace<Group2SparseVectorMatrix>>(5);
  TestRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolverInPlace<Group3SparseVectorMatrix>>(5);
  TestRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolverInPlace<Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolverInPlace, VectorOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto initial_value : initial_values)
  {
    TestExtremeInitialValue<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::LinearSolverInPlace<Group1SparseVectorMatrix>>(1, initial_value);
    TestExtremeInitialValue<
        Group2VectorMatrix,
        Group2SparseVectorMatrix,
        micm::LinearSolverInPlace<Group2SparseVectorMatrix>>(2, initial_value);
    TestExtremeInitialValue<
        Group3VectorMatrix,
        Group3SparseVectorMatrix,
        micm::LinearSolverInPlace<Group3SparseVectorMatrix>>(5, initial_value);
    TestExtremeInitialValue<
        Group4VectorMatrix,
        Group4SparseVectorMatrix,
        micm::LinearSolverInPlace<Group4SparseVectorMatrix>>(5, initial_value);
  }
}

TEST(LinearSolverInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolverInPlace<Group1SparseVectorMatrix>>(5);
  TestDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolverInPlace<Group2SparseVectorMatrix>>(5);
  TestDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolverInPlace<Group3SparseVectorMatrix>>(5);
  TestDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolverInPlace<Group4SparseVectorMatrix>>(5);
}

TEST(LinearSolverInPlace, VectorDiagonalMarkowitzReordering)
{
  TestMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  TestMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  TestMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  TestMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
