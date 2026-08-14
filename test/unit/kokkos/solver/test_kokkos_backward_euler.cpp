#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/process_set.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/kokkos/solver/kokkos_solver_builder.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

namespace
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::StoichSpecies(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::StoichSpecies(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 })
                         .SetPhase(gas_phase)
                         .Build();

  auto the_system = micm::System(gas_phase);
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(BackwardEuler, CanCallSolve)
{
  auto params = micm::BackwardEulerSolverParameters();

  auto be = micm::KokkosSolverBuilder<micm::BackwardEulerSolverParameters>(params)
                .SetSystem(the_system)
                .SetReactions(reactions)
                .Build();
  micm::Real time_step = 1.0;

  auto state = be.GetState(1);
  state.SetAbsoluteTolerances({ 1e-6, 1e-6, 1e-6 });

  state.variables_[0] = { 1.0, 0.0, 0.0 };
  state.conditions_[0].temperature_ = 272.5;
  state.conditions_[0].pressure_ = 101253.3;
  state.conditions_[0].air_density_ = 1e6;
  be.UpdateStateParameters(state);

  EXPECT_NO_THROW(auto result = be.Solve(time_step, state));
}

template<class DenseMatrixPolicy>
void CheckIsConverged()
{
  using LinearSolverPolicy = micm::LinearSolver<DenseMatrixPolicy, micm::StandardSparseMatrix>;
  using RatesPolicy = micm::ProcessSet<DenseMatrixPolicy, micm::StandardSparseMatrix>;
  using ConstraintSetPolicy = micm::ConstraintSet<DenseMatrixPolicy, micm::StandardSparseMatrix>;
  using BackwardEuler = micm::AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>;

  micm::BackwardEulerSolverParameters parameters;
  DenseMatrixPolicy residual{ 4, 3, 0.0 };
  DenseMatrixPolicy Yn1{ 4, 3, 0.0 };
  typename DenseMatrixPolicy::template ScalarType<micm::Bool> is_converged;

  parameters.small_ = 1e-6;
  micm::Real relative_tolerance = 1e-3;
  typename DenseMatrixPolicy::template VectorType<micm::Real> absolute_tolerance_data = { 1e-6, 1e-6, 1e-6 };

  residual.CopyToDevice();
  Yn1.CopyToDevice();
  absolute_tolerance_data.CopyToDevice();
  auto absolute_tolerance = std::as_const(absolute_tolerance_data).GetView();

  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_TRUE(is_converged);
  residual[0][1] = 1e-5;
  residual.CopyToDevice();
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_FALSE(is_converged);
  parameters.small_ = 1e-4;
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_TRUE(is_converged);
  residual[3][2] = 1e-3;
  residual.CopyToDevice();
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_FALSE(is_converged);
  Yn1[3][2] = 10.0;
  Yn1.CopyToDevice();
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_TRUE(is_converged);
  residual[3][2] = 1e-1;
  residual.CopyToDevice();
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_FALSE(is_converged);
  absolute_tolerance_data[2] = 1.0;
  absolute_tolerance_data.CopyToDevice();
  BackwardEuler::IsConverged(parameters, residual, Yn1, absolute_tolerance, relative_tolerance, is_converged);
  ASSERT_TRUE(is_converged);
}

TEST(BackwardEuler, IsConverged)
{
  CheckIsConverged<micm::KokkosDenseMatrix<micm::Real, 1>>();
  CheckIsConverged<micm::KokkosDenseMatrix<micm::Real, 2>>();
  CheckIsConverged<micm::KokkosDenseMatrix<micm::Real, 3>>();
  CheckIsConverged<micm::KokkosDenseMatrix<micm::Real, 4>>();
  CheckIsConverged<micm::KokkosDenseMatrix<micm::Real, 5>>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}