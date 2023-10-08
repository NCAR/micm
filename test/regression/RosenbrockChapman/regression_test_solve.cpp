#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>

#include "regression_test_solve_policy.hpp"

TEST(RegressionRosenbrock, TwoStageSolve)
{
  auto solver = getTwoStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testSolve(solver, 1.0e-2);
}

TEST(RegressionRosenbrock, ThreeStageSolve)
{
  auto solver = getThreeStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testSolve(solver);
}

TEST(RegressionRosenbrock, FourStageSolve)
{
  auto solver = getFourStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, FourStageDASolve)
{
  auto solver = getFourStageDAMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, SixStageDASolve)
{
  auto solver = getSixStageDAMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, VectorSolve)
{
  auto solver1 = getThreeStageMultiCellChapmanSolver<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::LinearSolver<double, Group1SparseVectorMatrix>>(3);
  testSolve(solver1);

  auto solver2 = getThreeStageMultiCellChapmanSolver<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::LinearSolver<double, Group2SparseVectorMatrix>>(3);
  testSolve(solver2);

  auto solver3 = getThreeStageMultiCellChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::LinearSolver<double, Group3SparseVectorMatrix>>(3);
  testSolve(solver3);

  auto solver4 = getThreeStageMultiCellChapmanSolver<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::LinearSolver<double, Group4SparseVectorMatrix>>(3);
  testSolve(solver4);
}