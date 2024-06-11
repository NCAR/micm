#include "regression_test_dforce_dy_policy.hpp"

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<std::size_t L>
using VectorBuilder = micm::CpuSolverBuilder<
      micm::RosenbrockSolverParameters,
      micm::VectorMatrix<double, L>,
      micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RegressionRosenbrock, Jacobian)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  auto solver = getChapmanSolver(builder, 3);
  testJacobian(solver);
}

TEST(RegressionRosenbrock, VectorJacobian)
{
  {
    auto builder = VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder, 3);
    testJacobian(solver);
  }
  {
    auto builder = VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder, 3);
    testJacobian(solver);
  }
  {
    auto builder = VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder, 3);
    testJacobian(solver);
  }
  {
    auto builder = VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
    auto solver = getChapmanSolver(builder, 3);
    testJacobian(solver);
  }
}