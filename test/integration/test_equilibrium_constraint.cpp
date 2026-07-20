// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace micm;

/// @brief Helper function to compute temperature-dependent equilibrium constant using Van't Hoff equation
micm::Real ComputeEquilibriumConstant(micm::Real K_HLC_ref, micm::Real delta_H, micm::Real T)
{
  return K_HLC_ref * std::exp((delta_H / constants::GAS_CONSTANT) * (1.0 / T - 1.0 / 298.15));
}

/// @brief This test verifies that the SolverBuilder correctly accepts constraints via SetConstraints
///        and that the resulting solver can be built without errors.
TEST(EquilibriumIntegration, SetConstraintsAPIWorks)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  micm::Real k = 0.1;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // C is an algebraic variable (not in any kinetic reaction)
  micm::Real K_eq = 0.034;
  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = -2400.0 }));

  // Build solver with constraints
  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify constraint metadata
  ASSERT_EQ(state.state_size_, 3);  // A, B, and C
  ASSERT_EQ(state.constraint_size_, 1);
  ASSERT_TRUE(state.variable_map_.contains("A"));
  ASSERT_TRUE(state.variable_map_.contains("B"));
  ASSERT_TRUE(state.variable_map_.contains("C"));
  ASSERT_TRUE(state.custom_rate_parameter_map_.contains("B_C_eq"));

  // Verify mass-matrix diagonal.
  // Size is species only; constrained species rows are algebraic (0), others are ODE (1).
  ASSERT_EQ(state.upper_left_identity_diagonal_.size(), 3);
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("A")], 1.0);  // A is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("B")], 1.0);  // B is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("C")], 0.0);  // C is algebraic

  // Verify state can be initialized
  micm::Index A_idx = state.variable_map_.at("A");
  micm::Index B_idx = state.variable_map_.at("B");
  micm::Index C_idx = state.variable_map_.at("C");
  micm::Index B_C_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.1;
  state.variables_[0][C_idx] = K_eq * 0.1;  // C = K_eq * B
  state.conditions_[0].temperature_ = 300.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters works
  solver.UpdateStateParameters(state);

  micm::Real expected_K_eq = ComputeEquilibriumConstant(K_eq, -2400.0, 300.0);
  // Scale the tolerance by machine epsilon so double keeps ~1e-16 strictness while float passes.
  EXPECT_NEAR(
      state.custom_rate_parameters_[0][B_C_eq_idx],
      expected_K_eq,
      std::abs(expected_K_eq) * 1e2 * std::numeric_limits<micm::Real>::epsilon());
}

/// @brief Verifies that multiple constraints can be added via SetConstraints and the solver
TEST(EquilibriumIntegration, SetConstraintsAPIMultipleConstraints)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");
  auto D = Species("D");
  auto E = Species("E");
  auto F = Species("F");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C, D, E, F } };

  // Simple kinetic reactions
  Process rxn1 = ChemicalReactionBuilder()
                     .SetReactants({ A })
                     .SetProducts({ { B, 1 } })
                     .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 0.5, .B_ = 0, .C_ = 0 })
                     .SetPhase(gas_phase)
                     .Build();

  Process rxn2 = ChemicalReactionBuilder()
                     .SetReactants({ D })
                     .SetProducts({ { E, 1 } })
                     .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 0.2, .B_ = 0, .C_ = 0 })
                     .SetPhase(gas_phase)
                     .Build();

  // Two equilibrium constraints
  micm::Real K_eq1 = 0.034;  // CO2-like equilibrium
  micm::Real K_eq2 = 0.012;  // Different gas equilibrium
  micm::Real delta_H1 = -2400.0;
  micm::Real delta_H2 = -2000.0;

  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq1, .delta_H_ = delta_H1 }));
  constraints.emplace_back(EquilibriumConstraint(
      "E_F_eq",
      F,
      std::vector<StoichSpecies>{ { E, 1.0 } },
      std::vector<StoichSpecies>{ { F, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq2, .delta_H_ = delta_H2 }));

  // Build solver with multiple constraints
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn1, rxn2 })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  ASSERT_EQ(state.state_size_, 6);       // A, B, C, D, E, F
  ASSERT_EQ(state.constraint_size_, 2);  // Two constraints

  ASSERT_EQ(state.upper_left_identity_diagonal_.size(), 6);
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("A")], 1.0);  // A is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("B")], 1.0);  // B is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("C")], 0.0);  // C is algebraic
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("D")], 1.0);  // D is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("E")], 1.0);  // E is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("F")], 0.0);  // F is algebraic

  // Initialize state and verify UpdateStateParameters works
  state.variables_[0][state.variable_map_.at("A")] = 1.0;
  state.variables_[0][state.variable_map_.at("B")] = 0.1;
  state.variables_[0][state.variable_map_.at("C")] = K_eq1 * 0.1;  // C should satisfy C = K_eq1 * B
  state.variables_[0][state.variable_map_.at("D")] = 0.5;
  state.variables_[0][state.variable_map_.at("E")] = 0.05;
  state.variables_[0][state.variable_map_.at("F")] = K_eq2 * 0.05;  // F should satisfy F = K_eq2 * E

  micm::Real current_temp = 310.0;
  state.conditions_[0].temperature_ = current_temp;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  // Verify temperature-dependent K_eq values are calculated correctly
  micm::Real expected_K_eq1 = ComputeEquilibriumConstant(K_eq1, delta_H1, current_temp);
  micm::Real expected_K_eq2 = ComputeEquilibriumConstant(K_eq2, delta_H2, current_temp);
  EXPECT_NEAR(state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")], expected_K_eq1, 1e-10);
  // Scale the tolerance by machine epsilon so double keeps ~1e-16 strictness while float passes.
  EXPECT_NEAR(
      state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("E_F_eq")],
      expected_K_eq2,
      std::abs(expected_K_eq2) * 1e2 * std::numeric_limits<micm::Real>::epsilon());
}

/// @brief Test DAE solving - actually calls Solve() with algebraic constraints
/// This exercises mass-matrix DAE enforcement where constrained species rows are algebraic.
TEST(EquilibriumIntegration, DAESolveWithConstraint)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  micm::Real k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  // This couples B (ODE variable) to C (algebraic variable)
  micm::Real K_eq = 2.0;
  micm::Real delta_H = -2400.0;
  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = delta_H }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Float precision cannot advance the default initial internal step (DEFAULT_H_START * time_step)
  // once simulated time is O(1) (the step falls below the float ULP); start with an absolute step
  // above unit round-off. Double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
    options.h_start_ = 1.0e-6;
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  micm::Index A_idx = state.variable_map_.at("A");
  micm::Index B_idx = state.variable_map_.at("B");
  micm::Index C_idx = state.variable_map_.at("C");
  micm::Index B_C_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  // Initial conditions: A=1, B=0, C=0
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 270.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly before time integration
  solver.UpdateStateParameters(state);
  micm::Real expected_K_eq = ComputeEquilibriumConstant(K_eq, delta_H, 270.0);
  // Scale the tolerance by machine epsilon so double keeps ~1e-16 strictness while float passes.
  EXPECT_NEAR(
      state.custom_rate_parameters_[0][B_C_eq_idx],
      expected_K_eq,
      std::abs(expected_K_eq) * 1e2 * std::numeric_limits<micm::Real>::epsilon());

  // Solve with smaller time steps
  micm::Real dt = 0.001;
  micm::Real total_time = 0.1;
  micm::Real time = 0.0;
  micm::Index steps = 0;

  while (time < total_time)
  {
    // Updates temperature-dependent state parameters.
    // Because state conditions remains constant within the while loop,
    // these values do not change during execution.
    // This behavior may be revised in the future.
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);

    if (result.state_ != SolverState::Converged)
    {
      FAIL() << "DAE solve did not converge at step " << steps << ", time=" << time;
    }

    // Verify constraint is maintained by the solver
    micm::Real constraint_residual =
        state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(constraint_residual, 0.0, 1.0e-6)
        << "Constraint not satisfied at step " << steps << ": K_eq*B - C = " << constraint_residual;

    time += dt;
    steps++;
  }

  // Verify constraint is satisfied: C = K_eq * B
  micm::Real expected_C = state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx];
  micm::Real final_residual =
      state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];

  EXPECT_NEAR(state.variables_[0][C_idx], expected_C, 1.0e-6);
  EXPECT_NEAR(final_residual, 0.0, 1.0e-6);

  // Verify mass conservation: A + B should be conserved (approximately)
  micm::Real total = state.variables_[0][A_idx] + state.variables_[0][B_idx];
  EXPECT_NEAR(total, 1.0, 0.01);
}

/// @brief Test DAE solve with constraints and state reordering enabled
TEST(EquilibriumIntegration, DAESolveWithConstraintAndReorderState)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  micm::Real k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  micm::Real K_eq = 2.0;
  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Float precision cannot advance the default initial internal step (DEFAULT_H_START * time_step)
  // once simulated time is O(1) (the step falls below the float ULP); start with an absolute step
  // above unit round-off. Double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
    options.h_start_ = 1.0e-6;
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(true)
                    .Build();

  auto state = solver.GetState(1);

  micm::Index A_idx = state.variable_map_.at("A");
  micm::Index B_idx = state.variable_map_.at("B");
  micm::Index C_idx = state.variable_map_.at("C");
  micm::Index B_C_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 400.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly
  solver.UpdateStateParameters(state);
  micm::Real expected_K_eq = ComputeEquilibriumConstant(K_eq, -2400.0, 400.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")], expected_K_eq, 1e-10);

  micm::Real dt = 0.001;
  micm::Real total_time = 0.1;
  micm::Real time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Reordered DAE solve did not converge at time=" << time;

    // Constraint should hold at each step
    micm::Real residual = state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-6) << "Constraint violated at time=" << time;

    time += dt;
  }

  // Mass conservation for A -> B
  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test DAE solve with two coupled constraints sharing a species
/// A -> B (kinetic), B <-> C (constraint 1), B <-> D (constraint 2)
TEST(EquilibriumIntegration, DAESolveWithTwoCoupledConstraints)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");
  auto D = Species("D");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C, D } };

  micm::Real k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  micm::Real K_eq1 = 3.0;
  micm::Real K_eq2 = 5.0;
  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq1, .delta_H_ = -2400.0 }));
  constraints.emplace_back(EquilibriumConstraint(
      "B_D_eq",
      D,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { D, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq2, .delta_H_ = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Float precision cannot advance the default initial internal step (DEFAULT_H_START * time_step)
  // once simulated time is O(1) (the step falls below the float ULP); start with an absolute step
  // above unit round-off. Double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
    options.h_start_ = 1.0e-6;
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  micm::Index A_idx = state.variable_map_.at("A");
  micm::Index B_idx = state.variable_map_.at("B");
  micm::Index C_idx = state.variable_map_.at("C");
  micm::Index D_idx = state.variable_map_.at("D");
  micm::Index B_C_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");
  micm::Index B_D_eq_idx = state.custom_rate_parameter_map_.at("B_D_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.variables_[0][D_idx] = 0.0;
  state.conditions_[0].temperature_ = 300.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly for both constraints
  solver.UpdateStateParameters(state);
  micm::Real expected_K_eq1 = ComputeEquilibriumConstant(K_eq1, -2400.0, 300.0);
  micm::Real expected_K_eq2 = ComputeEquilibriumConstant(K_eq2, -2400.0, 300.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_C_eq_idx], expected_K_eq1, 1e-10);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_D_eq_idx], expected_K_eq2, 1e-10);

  micm::Real dt = 0.001;
  micm::Real total_time = 0.1;
  micm::Real time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Coupled constraints did not converge at time=" << time;

    micm::Real residual1 =
        state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    micm::Real residual2 =
        state.custom_rate_parameters_[0][B_D_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][D_idx];
    EXPECT_NEAR(residual1, 0.0, 1.0e-6) << "Constraint 1 violated at time=" << time;
    EXPECT_NEAR(residual2, 0.0, 1.0e-6) << "Constraint 2 violated at time=" << time;

    time += dt;
  }

  // Both constraints should be satisfied at the end
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx], state.variables_[0][C_idx], 1.0e-6);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_D_eq_idx] * state.variables_[0][B_idx], state.variables_[0][D_idx], 1.0e-6);
  // Mass conservation: A + B should be conserved
  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test DAE solve with non-unit stoichiometric coefficient
/// 2A <-> B means K_eq * [A]^2 - [B] = 0
TEST(EquilibriumIntegration, DAESolveWithNonUnitStoichiometry)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Reaction: C -> A (to produce A for the equilibrium)
  micm::Real k = 0.5;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ C })
                    .SetProducts({ { A, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * [A]^2 - [B] = 0
  micm::Real K_eq = 10.0;
  std::vector<Constraint> constraints;
  constraints.emplace_back(EquilibriumConstraint(
      "A2_B_eq",
      B,
      std::vector<StoichSpecies>{ { A, 2.0 } },
      std::vector<StoichSpecies>{ { B, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Float precision cannot advance the default initial internal step (DEFAULT_H_START * time_step)
  // once simulated time is O(1) (the step falls below the float ULP); start with an absolute step
  // above unit round-off. Double mode keeps the original default (h_start_ == 0.0).
  if constexpr (!std::is_same_v<micm::Real, double>)
    options.h_start_ = 1.0e-6;
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  micm::Index A_idx = state.variable_map_.at("A");
  micm::Index B_idx = state.variable_map_.at("B");
  micm::Index C_idx = state.variable_map_.at("C");
  micm::Index A2_B_eq_idx = state.custom_rate_parameter_map_.at("A2_B_eq");

  // Start with some A so the constraint has something to work with
  state.variables_[0][A_idx] = 0.1;
  state.variables_[0][B_idx] = K_eq * 0.1 * 0.1;  // B = K_eq * A^2
  state.variables_[0][C_idx] = 1.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly
  solver.UpdateStateParameters(state);
  micm::Real expected_K_eq = ComputeEquilibriumConstant(K_eq, -2400.0, 298.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][A2_B_eq_idx], expected_K_eq, 1e-10);

  micm::Real dt = 0.001;
  micm::Real total_time = 0.05;
  micm::Real time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "NonUnit stoich did not converge at time=" << time;

    // Constraint: K_eq * [A]^2 - [B] = 0
    micm::Real A_val = state.variables_[0][A_idx];
    micm::Real B_val = state.variables_[0][B_idx];
    micm::Real residual = state.custom_rate_parameters_[0][A2_B_eq_idx] * A_val * A_val - B_val;
    EXPECT_NEAR(residual, 0.0, 1.0e-5) << "Constraint violated at time=" << time;

    time += dt;
  }
}