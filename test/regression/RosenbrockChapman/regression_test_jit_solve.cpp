#include "jit_util.hpp"
#include "regression_test_solve_policy.hpp"

#include <micm/solver/jit_linear_solver.hpp>
#include <micm/solver/jit_rosenbrock.hpp>

#include <gtest/gtest.h>

TEST(RegressionJitRosenbrock, TwoStageSolve)
{
  auto solver = getTwoStageMultiCellJitChapmanSolver<3>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, ThreeStageSolve)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<3>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, FourStageSolve)
{
  auto solver = getFourStageMultiCellJitChapmanSolver<3>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, FourStageDASolve)
{
  auto solver = getFourStageDAMultiCellJitChapmanSolver<3>(3);
  testSolve<>(solver, 1.0e-2);
}

TEST(RegressionJitRosenbrock, SixStageDASolve)
{
  auto solver = getSixStageDAMultiCellJitChapmanSolver<3>(3);
  testSolve<>(solver, 1.0e-2);
}