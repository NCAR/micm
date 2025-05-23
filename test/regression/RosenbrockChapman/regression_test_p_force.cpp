#include "regression_test_p_force_policy.hpp"

#include <gtest/gtest.h>

template<std::size_t L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, RateConstants)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testRateConstants(solver);
}

TEST(RegressionRosenbrock, VectorRateConstants)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testRateConstants(solver);
  }
}

TEST(RegressionRosenbrock, Forcing)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testForcing<micm::Matrix<double>>(solver);
}

TEST(RegressionRosenbrock, VectorForcing)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testForcing<micm::VectorMatrix<double, 1>>(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testForcing<micm::VectorMatrix<double, 2>>(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testForcing<micm::VectorMatrix<double, 3>>(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testForcing<micm::VectorMatrix<double, 4>>(solver);
  }
}