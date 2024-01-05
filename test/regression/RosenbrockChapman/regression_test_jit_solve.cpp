#include <gtest/gtest.h>

#include <micm/solver/jit_linear_solver.hpp>
#include <micm/solver/jit_rosenbrock.hpp>

#include "jit_util.hpp"
#include "regression_test_solve_policy.hpp"

TEST(RegressionJitRosenbrock, TwoStageSolve)
{
  auto solver = getTwoStageMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, ThreeStageSolve)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, FourStageSolve)
{
  auto solver = getFourStageMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, FourStageDASolve)
{
  auto solver = getFourStageDAMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, SixStageDASolve)
{
  auto solver = getSixStageDAMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testSolve<>(solver, 1.0e-2);
}