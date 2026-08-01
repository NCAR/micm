#include "regression_test_p_force_policy.hpp"

#include <micm/util/types.hpp>

#include <gtest/gtest.h>

template<micm::Index L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, RateConstants)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestRateConstants(solver);
}

TEST(RegressionRosenbrock, VectorRateConstants)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestRateConstants(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestRateConstants(solver);
  }
}

TEST(RegressionRosenbrock, Forcing)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestForcing<micm::Matrix<micm::Real>>(solver);
}

TEST(RegressionRosenbrock, VectorForcing)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestForcing<micm::VectorMatrix<micm::Real, 1>>(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestForcing<micm::VectorMatrix<micm::Real, 2>>(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestForcing<micm::VectorMatrix<micm::Real, 3>>(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestForcing<micm::VectorMatrix<micm::Real, 4>>(solver);
  }
}