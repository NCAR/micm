#include <micm/process/process_set.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

namespace
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yields(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yields(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase);

  auto the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(BackwardEuler, CanCallSolve)
{
  auto params = micm::BackwardEulerSolverParameters();
  params.absolute_tolerance_ = { 1e-6, 1e-6, 1e-6 };

  auto be = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(params)
                .SetSystem(the_system)
                .SetReactions(reactions)
                .SetNumberOfGridCells(1)
                .Build();
  double time_step = 1.0;

  auto state = be.GetState();

  state.variables_[0] = { 1.0, 0.0, 0.0 };
  state.conditions_[0].temperature_ = 272.5;
  state.conditions_[0].pressure_ = 101253.3;
  state.conditions_[0].air_density_ = 1e6;
  be.CalculateRateConstants(state);

  EXPECT_NO_THROW(auto result = be.Solve(time_step, state));
}

template<class DenseMatrixPolicy>
void CheckIsConverged()
{
  using LinearSolverPolicy = micm::LinearSolver<micm::StandardSparseMatrix>;
  using RatesPolicy = micm::ProcessSet;
  using BackwardEuler = micm::BackwardEuler<RatesPolicy, LinearSolverPolicy>;

  micm::BackwardEulerSolverParameters parameters;
  DenseMatrixPolicy residual{ 4, 3, 0.0 };
  DenseMatrixPolicy state{ 4, 3, 0.0 };

  parameters.small_ = 1e-6;
  parameters.relative_tolerance_ = 1e-3;
  parameters.absolute_tolerance_ = { 1e-6, 1e-6, 1e-6 };

  ASSERT_TRUE(BackwardEuler::IsConverged(parameters, residual, state));
  residual[0][1] = 1e-5;
  ASSERT_FALSE(BackwardEuler::IsConverged(parameters, residual, state));
  parameters.small_ = 1e-4;
  ASSERT_TRUE(BackwardEuler::IsConverged(parameters, residual, state));
  residual[3][2] = 1e-3;
  ASSERT_FALSE(BackwardEuler::IsConverged(parameters, residual, state));
  state[3][2] = 10.0;
  ASSERT_TRUE(BackwardEuler::IsConverged(parameters, residual, state));
  residual[3][2] = 1e-1;
  ASSERT_FALSE(BackwardEuler::IsConverged(parameters, residual, state));
  parameters.absolute_tolerance_[2] = 1.0;
  ASSERT_TRUE(BackwardEuler::IsConverged(parameters, residual, state));
}

TEST(BackwardEuler, IsConverged)
{
  CheckIsConverged<micm::Matrix<double>>();
  CheckIsConverged<micm::VectorMatrix<double, 1>>();
  CheckIsConverged<micm::VectorMatrix<double, 2>>();
  CheckIsConverged<micm::VectorMatrix<double, 3>>();
  CheckIsConverged<micm::VectorMatrix<double, 4>>();
  CheckIsConverged<micm::VectorMatrix<double, 5>>();
}
