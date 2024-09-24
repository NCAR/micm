#include "../../solver/test_linear_solver_policy.hpp"

#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <functional>

using Group1VectorMatrix = micm::VectorMatrix<double, 1>;
using Group2VectorMatrix = micm::VectorMatrix<double, 2>;
using Group3VectorMatrix = micm::VectorMatrix<double, 3>;
using Group4VectorMatrix = micm::VectorMatrix<double, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(JitLinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::JitLinearSolver<1, Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>>();
}

TEST(JitLinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::JitLinearSolver<1, Group1SparseVectorMatrix>>(1);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitLinearSolver<2, Group2SparseVectorMatrix>>(2);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::JitLinearSolver<3, Group3SparseVectorMatrix>>(3);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::JitLinearSolver<4, Group4SparseVectorMatrix>>(4);
}

TEST(JitLinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::JitLinearSolver<1, Group1SparseVectorMatrix>>(1);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitLinearSolver<2, Group2SparseVectorMatrix>>(2);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::JitLinearSolver<3, Group3SparseVectorMatrix>>(3);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::JitLinearSolver<4, Group4SparseVectorMatrix>>(4);
}

TEST(JitLinearSolver, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, INFINITY };
  for (auto initial_value : initial_values)
  {
    testExtremeInitialValue<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::JitLinearSolver<1, Group1SparseVectorMatrix>>(1, initial_value);
  }
}