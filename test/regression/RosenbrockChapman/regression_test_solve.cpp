#include "regression_test_solve_policy.hpp"

#include <micm/util/types.hpp>

#include <gtest/gtest.h>
template<micm::Index L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, TwoStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestSolve(solver, 1.0e-2);
}

TEST(RegressionRosenbrock, ThreeStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, FourStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, FourStageDASolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, SixStageDASolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, VectorSolve)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestSolve(solver, 1.0e-4);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestSolve(solver, 1.0e-4);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestSolve(solver, 1.0e-4);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestSolve(solver, 1.0e-4);
  }
}
