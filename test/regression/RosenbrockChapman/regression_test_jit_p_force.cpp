#include "jit_util.hpp"
#include "regression_test_p_force_policy.hpp"

#include <micm/solver/jit_linear_solver.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

TEST(RegressionJitRosenbrock, VectorRateConstants)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<3>(3);
  testRateConstants<>(solver);
}

TEST(RegressionJitRosenbrock, VectorForcing)
{
  auto solver = getThreeStageMultiCellJitChapmanSolver<3>(3);
  testForcing<micm::VectorMatrix<double, 3>>(solver);
}