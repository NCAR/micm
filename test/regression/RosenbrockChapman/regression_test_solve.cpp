#include "regression_test_solve_policy.hpp"

#include <gtest/gtest.h>
template<std::size_t L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, TwoStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testSolve(solver, 1.0e-2);
}

TEST(RegressionRosenbrock, ThreeStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testSolve(solver);
}

TEST(RegressionRosenbrock, FourStageSolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::FourStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, FourStageDASolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, SixStageDASolve)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  auto solver = getChapmanSolver(builder);
  testSolve(solver, 1.0e-4);
}

TEST(RegressionRosenbrock, VectorSolve)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testSolve(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testSolve(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testSolve(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder);
    testSolve(solver);
  }
}