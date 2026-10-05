#include "regression_test_dforce_dy_policy.hpp"

#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

template<micm::Index L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, Jacobian)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = GetChapmanSolver(builder);
  TestJacobian(solver);
}

TEST(RegressionRosenbrock, VectorJacobian)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestJacobian(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestJacobian(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestJacobian(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = GetChapmanSolver(builder);
    TestJacobian(solver);
  }
}