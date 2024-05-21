#include "regression_test_p_force_policy.hpp"

#include <micm/solver/rosenbrock.hpp>

#include <gtest/gtest.h>

TEST(RegressionRosenbrock, RateConstants)
{
  auto solver = getThreeStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<SparseMatrix>>(3);
  testRateConstants<>(solver);
}

TEST(RegressionRosenbrock, VectorRateConstants)
{
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::LinearSolver<Group1SparseVectorMatrix>>(3);
    testRateConstants<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group2VectorMatrix,
        Group2SparseVectorMatrix,
        micm::LinearSolver<Group2SparseVectorMatrix>>(3);
    testRateConstants<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group3VectorMatrix,
        Group3SparseVectorMatrix,
        micm::LinearSolver<Group3SparseVectorMatrix>>(3);
    testRateConstants<>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group4VectorMatrix,
        Group4SparseVectorMatrix,
        micm::LinearSolver<Group4SparseVectorMatrix>>(3);
    testRateConstants<>(solver);
  }
}

TEST(RegressionRosenbrock, Forcing)
{
  auto solver = getThreeStageMultiCellChapmanSolver<DenseMatrix, SparseMatrix, micm::LinearSolver<SparseMatrix>>(3);
  testForcing<DenseMatrix>(solver);
}

TEST(RegressionRosenbrock, VectorForcing)
{
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::LinearSolver<Group1SparseVectorMatrix>>(3);
    testForcing<Group1VectorMatrix>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group2VectorMatrix,
        Group2SparseVectorMatrix,
        micm::LinearSolver<Group2SparseVectorMatrix>>(3);
    testForcing<Group2VectorMatrix>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group3VectorMatrix,
        Group3SparseVectorMatrix,
        micm::LinearSolver<Group3SparseVectorMatrix>>(3);
    testForcing<Group3VectorMatrix>(solver);
  }
  {
    auto solver = getThreeStageMultiCellChapmanSolver<
        Group4VectorMatrix,
        Group4SparseVectorMatrix,
        micm::LinearSolver<Group4SparseVectorMatrix>>(3);
    testForcing<Group4VectorMatrix>(solver);
  }
}