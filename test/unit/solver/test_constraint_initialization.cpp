// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

using namespace micm;

/// @brief Helper: build a simple A→B system with equilibrium constraint C = K_eq * B
struct SimpleConstrainedSystem
{
  static constexpr micm::Real K = 0.5;
  static constexpr micm::Real K_EQ = 2.0;
  static constexpr micm::Real DELTA_H = 0.0;  // No temperature dependence for simplicity

  template<class SolverBuilderPolicy>
  static auto Build(SolverBuilderPolicy builder)
  {
    using DenseMatrix = SolverBuilderPolicy::DenseMatrixPolicyType;
    using SparseMatrix = SolverBuilderPolicy::SparseMatrixPolicyType;
    auto A = Species("A");
    auto B = Species("B");
    auto C = Species("C");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { B, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = K, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
    std::vector<Constraint<DenseMatrix, SparseMatrix>> constraints;
    constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrix>(
        "B_C_eq",
        C,
        std::vector<StoichSpecies>{ { B, 1.0 } },
        std::vector<StoichSpecies>{ { C, 1.0 } },
        { .K_HLC_ref_ = K_EQ, .delta_H_ = DELTA_H }));

    return builder.SetSystem(System(gas_phase))
        .SetReactions({ rxn })
        .SetConstraints(std::move(constraints))
        .SetReorderState(false)
        .Build();
  }
};

using StandardBuilder =
    CpuSolverBuilder<RosenbrockSolverParameters, Matrix<micm::Real>, SparseMatrix<micm::Real, SparseMatrixStandardOrdering>>;

/// @brief Projection-only system: Z is slaved to X by row_scale * (X - Z) = 0, so the
///        algebraic solution is Z = X for every row scale.
auto BuildScaledLinearConstraintSolver(const RosenbrockSolverParameters& parameters, micm::Real row_scale)
{
  using DenseMatrix = StandardBuilder::DenseMatrixPolicyType;
  using SparseMatrix = StandardBuilder::SparseMatrixPolicyType;

  const auto X = Species("X");
  const auto Z = Species("Z");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ X, Z } };

  std::vector<Constraint<DenseMatrix, SparseMatrix>> constraints;
  constraints.emplace_back(
      LinearConstraint<DenseMatrix, SparseMatrix>("scaled_copy", Z, { { X, row_scale }, { Z, -row_scale } }, 0.0));

  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .SetConstraints(std::move(constraints))
      .SetReorderState(false)
      .Build();
}

/// @brief Test that consistent initial conditions don't change state
TEST(ConstraintInitialization, ConsistentICsUnchanged)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_EQ * 0.5;  // C = K_eq * B = consistent
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);
  // A should not change due to initialization (only from time stepping)
  // B should not change due to initialization (only from time stepping)
  // C was already consistent, so initialization should be nearly a no-op
  // Constraint init should converge in ≤1 iteration
  EXPECT_LE(result.stats_.constraint_init_iterations_, 1);
}

/// @brief Test that mildly inconsistent ICs are corrected
TEST(ConstraintInitialization, MildlyInconsistentICsCorrected)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 5.0;  // Wrong: should be K_eq * B = 1.0
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  micm::Real A_before = state.variables_[0][A_idx];
  micm::Real B_before = state.variables_[0][B_idx];

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // ODE variable A should not have been modified by initialization
  // (it will change from time stepping, but the initialization should not touch it)
  // We can't check exact equality post-solve because time stepping changes A,
  // but we verify the initialization converged and the constraint is satisfied
  EXPECT_GT(result.stats_.constraint_init_iterations_, micm::Bool(false));

  // After solve, the constraint should be satisfied: C ≈ K_eq * B
  micm::Real K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  micm::Real residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Test severely inconsistent ICs converge
TEST(ConstraintInitialization, SeverelyInconsistentICsConverge)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 1000.0;  // Wildly off: should be K_eq * B = 1.0
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // Constraint should be satisfied after initialization + solve
  micm::Real K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  micm::Real residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Row-scale invariance must survive the public Solve() path, not just the projection.
///        A falsely converged projection hands the integrator a state off the manifold, which
///        it reports as StepSizeTooSmall rather than as a constraint failure.
TEST(ConstraintInitialization, ScaledConstraintRowsSolveIdenticallyThroughSolve)
{
  using DenseMatrix = StandardBuilder::DenseMatrixPolicyType;
  using SparseMatrix = StandardBuilder::SparseMatrixPolicyType;

  constexpr micm::Real TOTAL = 1.3;
  constexpr micm::Real A0 = 0.7;
  constexpr micm::Real B0 = 0.15;
  constexpr micm::Real C_CONSISTENT = TOTAL - A0 - B0;

  // A -> B with the conservation row A + B + C = TOTAL multiplied through by `scale`.
  // A + B is conserved by the reaction, so C stays at C_CONSISTENT once projected.
  auto solve_with_row_scale = [](micm::Real scale)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto C = Species("C");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { B, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 0.5, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    std::vector<Constraint<DenseMatrix, SparseMatrix>> constraints;
    constraints.emplace_back(LinearConstraint<DenseMatrix, SparseMatrix>(
        "mass", C, std::vector<StoichSpecies>{ { A, scale }, { B, scale }, { C, scale } }, scale * TOTAL));

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    if constexpr (!std::is_same_v<micm::Real, double>)
    {
      options.h_start_ = 1.0e-6;
    }
    auto solver = StandardBuilder(options)
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = A0;
    state.variables_[0][state.variable_map_.at("B")] = B0;
    state.variables_[0][state.variable_map_.at("C")] = 5.0;  // far off the manifold
    state.conditions_[0].temperature_ = 298.15;
    state.conditions_[0].pressure_ = 101325.0;

    solver.UpdateStateParameters(state);

    auto result = solver.Solve(0.001, state);
    return std::make_pair(result.state_, state.variables_[0][state.variable_map_.at("C")]);
  };

  const auto reference = solve_with_row_scale(1.0);
  EXPECT_EQ(reference.first, SolverState::Converged);
  EXPECT_NEAR(reference.second, C_CONSISTENT, 1.0e-5 * C_CONSISTENT);

  const micm::Real scales[] = { static_cast<micm::Real>(1.0e-12), static_cast<micm::Real>(1.0e12) };
  for (micm::Real scale : scales)
  {
    const auto scaled = solve_with_row_scale(scale);
    EXPECT_EQ(scaled.first, SolverState::Converged) << "constraint row scale " << scale;
    EXPECT_NEAR(scaled.second, reference.second, 1.0e-5 * C_CONSISTENT) << "constraint row scale " << scale;
  }
}

/// @brief Test pure ODE system (no constraints) is unaffected
TEST(ConstraintInitialization, PureODESystemUnaffected)
{
  auto A = Species("A");
  auto B = Species("B");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B } };

  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 0.5, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.variables_[0][state.variable_map_.at("A")] = 1.0;
  state.variables_[0][state.variable_map_.at("B")] = 0.0;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  auto result = solver.Solve(0.01, state);

  EXPECT_EQ(result.state_, SolverState::Converged);
  // No constraint initialization should have happened
  EXPECT_EQ(result.stats_.constraint_init_iterations_, micm::Bool(false));
}

/// @brief Test multi-cell systems with different inconsistency levels
TEST(ConstraintInitialization, MultiCellSystems)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(3);  // 3 grid cells

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  micm::Real K_eq = SimpleConstrainedSystem::K_EQ;

  // Cell 0: consistent
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = K_eq * 0.5;

  // Cell 1: mildly inconsistent
  state.variables_[1][A_idx] = 2.0;
  state.variables_[1][B_idx] = 0.3;
  state.variables_[1][C_idx] = 5.0;  // Should be K_eq * 0.3 = 0.6

  // Cell 2: severely inconsistent (but non-negative — equilibrium constraint clamps negative concentrations)
  state.variables_[2][A_idx] = 0.5;
  state.variables_[2][B_idx] = 0.1;
  state.variables_[2][C_idx] = 100.0;  // Should be K_eq * 0.1 = 0.2

  for (micm::Index i = 0; i < 3; ++i)
  {
    state.conditions_[i].temperature_ = 298.15;
    state.conditions_[i].pressure_ = 101325.0;
  }

  solver.UpdateStateParameters(state);

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // Verify constraint satisfied in all cells after solve
  micm::Real K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  for (micm::Index i = 0; i < 3; ++i)
  {
    micm::Real residual = K_eq_actual * state.variables_[i][B_idx] - state.variables_[i][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-5) << "Constraint not satisfied in cell " << i;
  }
}

/// @brief Test that subsequent Solve() calls re-check and re-initialize if needed
TEST(ConstraintInitialization, SubsequentSolveCallsReinitialize)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // In float precision the default initial internal step (DEFAULT_H_START * time_step)
  // falls below the unit round-off (float epsilon ~1.2e-7), tripping the
  // step-size-too-small guard on the first iteration. Start with a larger step in
  // float mode; double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
  {
    options.h_start_ = 1.0e-6;
  }
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  // Start consistent
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_EQ * 0.5;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);
  auto result1 = solver.Solve(0.001, state);
  EXPECT_EQ(result1.state_, SolverState::Converged);

  // Perturb C off-manifold between solve calls (simulating external event)
  state.variables_[0][C_idx] = 999.0;

  solver.UpdateStateParameters(state);
  auto result2 = solver.Solve(0.001, state);
  EXPECT_EQ(result2.state_, SolverState::Converged);

  // Constraint should be re-satisfied after second solve
  micm::Real K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  micm::Real residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Test SolverStateToString for the new enum value
TEST(ConstraintInitialization, SolverStateToStringNewValue)
{
  EXPECT_EQ(SolverStateToString(SolverState::ConstraintInitializationFailed), "Constraint Initialization Failed");
}

/// @brief Projection-only system nonlinear in the algebraic variable: 2*[B] - [C]^2 = 0,
///        so C = 1 when B = 0.5. One Newton update from a distant guess lands nowhere near it.
auto BuildSquaredConstraintSolver(const RosenbrockSolverParameters& parameters)
{
  using DenseMatrix = StandardBuilder::DenseMatrixPolicyType;
  using SparseMatrix = StandardBuilder::SparseMatrixPolicyType;

  const auto B = Species("B");
  const auto C = Species("C");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ B, C } };

  std::vector<Constraint<DenseMatrix, SparseMatrix>> constraints;
  constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrix>(
      "squared",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 2.0 } },
      { .K_HLC_ref_ = 2.0, .delta_H_ = 0.0 }));

  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .SetConstraints(std::move(constraints))
      .SetReorderState(false)
      .Build();
}

/// @brief A failed projection must leave the caller's state exactly as it was passed in.
///        Otherwise the caller is handed a half-applied Newton iterate that is neither their
///        input nor a solution, with nothing in the return value to say so.
TEST(ConstraintInitialization, FailedInitializationRestoresCallerState)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_iterations_ = 1;  // too few updates to reach C = 1 from C = 100

  auto solver = BuildSquaredConstraintSolver(parameters);
  auto state = solver.GetState(1);

  const auto B_idx = state.variable_map_.at("B");
  const auto C_idx = state.variable_map_.at("C");
  constexpr micm::Real B_in = 0.5;
  constexpr micm::Real C_in = 100.0;
  state.variables_[0][B_idx] = B_in;
  state.variables_[0][C_idx] = C_in;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);
  state.variables_.CopyToDevice();

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
  state.variables_.CopyToHost();

  ASSERT_EQ(status, SolverState::ConstraintInitializationFailed);
  EXPECT_EQ(state.variables_[0][C_idx], C_in);
  EXPECT_EQ(state.variables_[0][B_idx], B_in);
}

/// @brief A converged projection still returns the corrected state, not the snapshot
TEST(ConstraintInitialization, SuccessfulInitializationKeepsTheCorrection)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = BuildSquaredConstraintSolver(parameters);
  auto state = solver.GetState(1);

  const auto B_idx = state.variable_map_.at("B");
  const auto C_idx = state.variable_map_.at("C");
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 100.0;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);
  state.variables_.CopyToDevice();

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
  state.variables_.CopyToHost();

  ASSERT_EQ(status, SolverState::Converged);
  // Converged means the remaining correction is a tenth of the state-error scale, which with the
  // default absolute tolerance of 1e-3 leaves C within ~1e-4 of the manifold, not within roundoff.
  EXPECT_NEAR(state.variables_[0][C_idx], 1.0, 1.0e-3);
}

/// @brief The weighted Newton-correction measure is invariant to complete constraint-row scaling
TEST(ConstraintInitialization, WeightedCorrectionIsInvariantToConstraintRowScaling)
{
  const std::array<micm::Real, 3> row_scales{ 1.0e-12, 1.0, 1.0e12 };

  for (const micm::Real row_scale : row_scales)
  {
    const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
    auto solver = BuildScaledLinearConstraintSolver(parameters, row_scale);
    auto state = solver.GetState(1);

    constexpr micm::Real rtol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;
    constexpr micm::Real atol = std::is_same_v<micm::Real, double> ? 1.0e-12 : 1.0e-6;
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<micm::Real>(state.state_size_, atol));

    const auto X = state.variable_map_.at("X");
    const auto Z = state.variable_map_.at("Z");
    state.variables_[0][X] = 1.0;
    state.variables_[0][Z] = 4.0;
    state.variables_.CopyToDevice();

    SolverStats stats;
    const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
    state.variables_.CopyToHost();

    EXPECT_EQ(status, SolverState::Converged) << "row scale=" << row_scale;
    EXPECT_NEAR(state.variables_[0][Z], 1.0, atol) << "row scale=" << row_scale;
  }
}

/// @brief constraint_init_tolerance is a dimensionless fraction of the state-error scale.
///        The same off-manifold state is already converged under loose tolerances and needs a
///        Newton update under tight ones. constraint_init_iterations_ counts iterations entered,
///        so it is one greater than the number of updates applied.
TEST(ConstraintInitialization, WeightedCorrectionUsesStateTolerances)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  EXPECT_EQ(parameters.constraint_init_tolerance_, micm::Real(0.1));

  constexpr micm::Real deviation = std::is_same_v<micm::Real, double> ? 5.0e-8 : 5.0e-6;
  constexpr micm::Real loose_rtol = std::is_same_v<micm::Real, double> ? 1.0e-6 : 1.0e-4;
  constexpr micm::Real loose_atol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-5;
  constexpr micm::Real tight_rtol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-6;
  constexpr micm::Real tight_atol = std::is_same_v<micm::Real, double> ? 1.0e-12 : 1.0e-8;

  auto loose_solver = BuildScaledLinearConstraintSolver(parameters, 1.0);
  auto loose_state = loose_solver.GetState(1);
  loose_state.SetRelativeTolerance(loose_rtol);
  loose_state.SetAbsoluteTolerances(std::vector<micm::Real>(loose_state.state_size_, loose_atol));
  const auto loose_X = loose_state.variable_map_.at("X");
  const auto loose_Z = loose_state.variable_map_.at("Z");
  loose_state.variables_[0][loose_X] = 1.0;
  loose_state.variables_[0][loose_Z] = 1.0 + deviation;
  const micm::Real loose_initial_Z = loose_state.variables_[0][loose_Z];
  loose_state.variables_.CopyToDevice();

  SolverStats loose_stats;
  const auto loose_status = loose_solver.solver_.InitializeConstraints(loose_state, parameters, loose_stats);
  loose_state.variables_.CopyToHost();

  EXPECT_EQ(loose_status, SolverState::Converged);
  EXPECT_EQ(loose_state.variables_[0][loose_Z], loose_initial_Z);
  EXPECT_EQ(loose_stats.constraint_init_iterations_, 1);  // no update applied

  auto tight_solver = BuildScaledLinearConstraintSolver(parameters, 1.0);
  auto tight_state = tight_solver.GetState(1);
  tight_state.SetRelativeTolerance(tight_rtol);
  tight_state.SetAbsoluteTolerances(std::vector<micm::Real>(tight_state.state_size_, tight_atol));
  const auto tight_X = tight_state.variable_map_.at("X");
  const auto tight_Z = tight_state.variable_map_.at("Z");
  tight_state.variables_[0][tight_X] = 1.0;
  tight_state.variables_[0][tight_Z] = 1.0 + deviation;
  tight_state.variables_.CopyToDevice();

  SolverStats tight_stats;
  const auto tight_status = tight_solver.solver_.InitializeConstraints(tight_state, parameters, tight_stats);
  tight_state.variables_.CopyToHost();

  EXPECT_EQ(tight_status, SolverState::Converged);
  EXPECT_NEAR(tight_state.variables_[0][tight_Z], 1.0, tight_atol);
  EXPECT_EQ(tight_stats.constraint_init_iterations_, 2);  // one update applied
}

/// @brief The three line-search controls ship with the documented defaults, and 24 backtracks is a
///        step length of 2^-24, below which a correction cannot move a double-precision iterate.
TEST(ConstraintInitialization, LineSearchParameterDefaults)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  EXPECT_EQ(parameters.constraint_init_max_backtracks_, micm::Index(24));
  EXPECT_EQ(parameters.constraint_init_backtrack_factor_, micm::Real(0.5));
  EXPECT_EQ(parameters.constraint_init_sufficient_decrease_, micm::Real(1.0e-4));
}

// ---------------------------------------------------------------------------
// Damped-Newton (line-search) regression tests
// ---------------------------------------------------------------------------

/// @brief Constraint-only external model enforcing row_scale * ([Z]^2 - [X]) = 0, so the algebraic
///        solution is Z = sqrt(X) for every row scale. The residual is nonlinear in the solved
///        variable, so an undamped Newton step from a guess far below the root overshoots by
///        1/(2*Z_0), which is what the line search exists to control. Multiplying the whole row by
///        row_scale changes G and dG/dy together, so it must change no decision the projection makes.
class ScaledSquareRootConstraintModel
{
 public:
  ScaledSquareRootConstraintModel(std::string x, std::string z, micm::Real row_scale)
      : x_(std::move(x)),
        z_(std::move(z)),
        row_scale_(row_scale)
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { z_ };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { x_, z_ };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    const auto i_x = state_indices.at(x_);
    const auto i_z = state_indices.at(z_);
    return { { i_z, i_x }, { i_z, i_z } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
  ConstraintUpdateStateParametersFunction(const std::unordered_map<std::string, micm::Index>&) const
  {
    return [](const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, micm::Index>&,
      const std::unordered_map<std::string, micm::Index>& var) const
  {
    const auto i_x = var.at(x_);
    const auto i_z = var.at(z_);
    const micm::Real scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (micm::Index i = 0; i < state.NumRows(); ++i)
      {
        forcing[i][i_z] = scale * (state[i][i_z] * state[i][i_z] - state[i][i_x]);
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, micm::Index>&,
      const std::unordered_map<std::string, micm::Index>& var,
      const SparseMatrixPolicy&) const
  {
    const auto i_x = var.at(x_);
    const auto i_z = var.at(z_);
    const micm::Real scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (micm::Index i = 0; i < jac.NumberOfBlocks(); ++i)
      {
        jac[i][i_z][i_z] -= 2.0 * scale * state[i][i_z];
        jac[i][i_z][i_x] -= -scale;
      }
    };
  }

 private:
  std::string x_;
  std::string z_;
  micm::Real row_scale_;
};

/// @brief Projection-only system nonlinear in the algebraic variable, with the whole constraint row
///        multiplied by row_scale: row_scale * ([Z]^2 - [X]) = 0, so Z = sqrt(X).
auto BuildScaledSquareRootSolver(const RosenbrockSolverParameters& parameters, micm::Real row_scale)
{
  const auto X = Species("X");
  const auto Z = Species("Z");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ X, Z } };

  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .AddExternalModel(ScaledSquareRootConstraintModel("X", "Z", row_scale))
      .SetReorderState(false)
      .Build();
}

/// @brief Drive the projection directly and report status, the converged Z, and the iteration count.
struct ProjectionOutcome
{
  SolverState status_;
  micm::Real z_;
  micm::Index iterations_;
  micm::Index solves_;
};

ProjectionOutcome ProjectSquareRoot(const RosenbrockSolverParameters& parameters, micm::Real row_scale, micm::Real z_initial)
{
  auto solver = BuildScaledSquareRootSolver(parameters, row_scale);
  auto state = solver.GetState(1);

  constexpr micm::Real rtol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;
  constexpr micm::Real atol = std::is_same_v<micm::Real, double> ? 1.0e-12 : 1.0e-6;
  state.SetRelativeTolerance(rtol);
  state.SetAbsoluteTolerances(std::vector<micm::Real>(state.state_size_, atol));

  state.variables_[0][state.variable_map_.at("X")] = 1.0;
  state.variables_[0][state.variable_map_.at("Z")] = z_initial;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);
  state.variables_.CopyToDevice();

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
  state.variables_.CopyToHost();

  return { status, state.variables_[0][state.variable_map_.at("Z")], stats.constraint_init_iterations_, stats.solves_ };
}

/// @brief Disabling the line search restores the undamped iteration exactly, so the change is
///        opt-out and bisectable. Both directions are pinned: the cold start must FAIL with the
///        search off, and the caller's state must come back untouched.
TEST(ConstraintInitialization, DisablingBacktracksReproducesUndampedBehavior)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_backtracks_ = 0;

  const auto outcome = ProjectSquareRoot(parameters, 1.0, 1.0e-6);

  EXPECT_EQ(outcome.status_, SolverState::ConstraintInitializationFailed);
  EXPECT_EQ(outcome.z_, micm::Real(1.0e-6));  // exact restoration, not approximate
  // No line-search solves were taken: one solve per Newton update and nothing more.
  EXPECT_EQ(outcome.solves_, outcome.iterations_);
}
