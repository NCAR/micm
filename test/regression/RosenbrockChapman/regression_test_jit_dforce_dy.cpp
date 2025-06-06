#include "jit_util.hpp"
#include "regression_test_dforce_dy_policy.hpp"

#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/jit/solver/jit_rosenbrock.hpp>

#include <gtest/gtest.h>

TEST(RegressionJitRosenbrock, VectorJacobian)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<3>();
  testJacobian<>(solver);
}