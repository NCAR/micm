#include "regression_test_dforce_dy_policy.hpp"

#include <micm/solver/rosenbrock.hpp>

#include <gtest/gtest.h>

TEST(RegressionRosenbrock, Jacobian)
{
  auto solver = getThreeStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<SparseMatrix>>(3);
  testJacobian<>(solver);
}

TEST(RegressionRosenbrock, VectorJacobian)
{
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::LinearSolver<Group1SparseVectorMatrix>>(3);
    testJacobian<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group2VectorMatrix,
        Group2SparseVectorMatrix,
        micm::LinearSolver<Group2SparseVectorMatrix>>(3);
    testJacobian<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group3VectorMatrix,
        Group3SparseVectorMatrix,
        micm::LinearSolver<Group3SparseVectorMatrix>>(3);
    testJacobian<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group4VectorMatrix,
        Group4SparseVectorMatrix,
        micm::LinearSolver<Group4SparseVectorMatrix>>(3);
    testJacobian<>(solver);
  }
}