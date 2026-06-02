#include "regression_policy.hpp"

#include <gtest/gtest.h>

auto rosenbrock_2stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
    micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());

auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters());

TEST(RosenbrockRegression, FlowTube)
{
  TestFlowTube(rosenbrock_2stage, "expected_results/rosenbrock_wall_loss.csv");
}

TEST(BackwardEulerRegression, FlowTube)
{
  TestFlowTube(backward_euler, "expected_results/backward_euler_wall_loss.csv");
}