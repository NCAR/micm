#include <gtest/gtest.h>

#include <micm/solver/jit_linear_solver.hpp>
#include <micm/solver/jit_rosenbrock.hpp>

#include "jit_util.hpp"
#include "regression_test_dforce_dy_policy.hpp"

TEST(RegressionJitRosenbrock, VectorJacobian)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::JitLinearSolver<3, Group3SparseVectorMatrix>,
      micm::JitProcessSet<3>>(3);
  testJacobian<>(solver);
}