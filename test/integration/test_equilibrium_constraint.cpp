// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace micm;

/// @brief Helper function to compute temperature-dependent equilibrium constant using Van't Hoff equation
double ComputeEquilibriumConstant(double K_HLC_ref, double delta_H, double T)
{
  return K_HLC_ref * std::exp((delta_H / constants::GAS_CONSTANT) * (1.0 / T - 1.0 / 298.15));
}

/// @brief This test verifies that the SolverBuilder correctly accepts constraints via SetConstraints
///        and that the resulting solver can be built without errors.
TEST(EquilibriumIntegration, SetConstraintsAPIWorks)
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  double k = 0.1;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // C is an algebraic variable (not in any kinetic reaction)
  double k_eq = 0.034;
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 }));

  // Build solver with constraints
  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify constraint metadata
  ASSERT_EQ(state.state_size_, 3);  // A, B, and C
  ASSERT_EQ(state.constraint_size_, 1);
  ASSERT_TRUE(state.variable_map_.count("A") > 0);
  ASSERT_TRUE(state.variable_map_.count("B") > 0);
  ASSERT_TRUE(state.variable_map_.count("C") > 0);
  ASSERT_TRUE(state.custom_rate_parameter_map_.count("B_C_eq") > 0);

  // Verify mass-matrix diagonal.
  // Size is species only; constrained species rows are algebraic (0), others are ODE (1).
  ASSERT_EQ(state.upper_left_identity_diagonal_.size(), 3);
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("A")], 1.0);  // A is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("B")], 1.0);  // B is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("C")], 0.0);  // C is algebraic

  // Verify state can be initialized
  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");
  std::size_t b_c_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.1;
  state.variables_[0][C_idx] = k_eq * 0.1;  // C = K_eq * B
  state.conditions_[0].temperature_ = 300.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters works
  solver.UpdateStateParameters(state);

  double expected_k_eq = ComputeEquilibriumConstant(k_eq, -2400.0, 300.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_C_eq_idx], expected_k_eq, 1e-10);
}

/// @brief Verifies that multiple constraints can be added via SetConstraints and the solver
TEST(EquilibriumIntegration, SetConstraintsAPIMultipleConstraints)
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");
  auto d = Species("D");
  auto e = Species("E");
  auto f = Species("F");

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
  double k_eq1 = 0.034;  // CO2-like equilibrium
  double k_eq2 = 0.012;  // Different gas equilibrium
  double delta_h1 = -2400.0;
  double delta_h2 = -2000.0;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq1, .delta_H = delta_H1 }));
  constraints.push_back(EquilibriumConstraint(
      "E_F_eq",
      F,
      std::vector<StoichSpecies>{ { E, 1.0 } },
      std::vector<StoichSpecies>{ { F, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq2, .delta_H = delta_H2 }));

  // Build solver with multiple constraints
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
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
  state.variables_[0][state.variable_map_.at("C")] = k_eq1 * 0.1;  // C should satisfy C = K_eq1 * B
  state.variables_[0][state.variable_map_.at("D")] = 0.5;
  state.variables_[0][state.variable_map_.at("E")] = 0.05;
  state.variables_[0][state.variable_map_.at("F")] = k_eq2 * 0.05;  // F should satisfy F = K_eq2 * E

  double current_temp = 310.0;
  state.conditions_[0].temperature_ = current_temp;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  // Verify temperature-dependent K_eq values are calculated correctly
  double expected_k_eq1 = ComputeEquilibriumConstant(k_eq1, delta_h1, current_temp);
  double expected_k_eq2 = ComputeEquilibriumConstant(k_eq2, delta_h2, current_temp);
  EXPECT_NEAR(state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")], expected_k_eq1, 1e-10);
  EXPECT_NEAR(state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("E_F_eq")], expected_k_eq2, 1e-10);
}

/// @brief Test DAE solving - actually calls Solve() with algebraic constraints
/// This exercises mass-matrix DAE enforcement where constrained species rows are algebraic.
TEST(EquilibriumIntegration, DAESolveWithConstraint)
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  double k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  // This couples B (ODE variable) to C (algebraic variable)
  double k_eq = 2.0;
  double delta_h = -2400.0;
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = delta_H }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");
  std::size_t b_c_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  // Initial conditions: A=1, B=0, C=0
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 270.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly before time integration
  solver.UpdateStateParameters(state);
  double expected_k_eq = ComputeEquilibriumConstant(k_eq, delta_h, 270.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_C_eq_idx], expected_k_eq, 1e-10);

  // Solve with smaller time steps
  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;
  int steps = 0;

  while (time < total_time)
  {
    // Updates temperature-dependent state parameters.
    // Because state conditions remains constant within the while loop,
    // these values do not change during execution.
    // This behavior may be revised in the future.
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);

    if (result.state_ != SolverState::CONVERGED)
    {
      FAIL() << "DAE solve did not converge at step " << steps << ", time=" << time;
    }

    // Verify constraint is maintained by the solver
    double constraint_residual =
        state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(constraint_residual, 0.0, 1.0e-6)
        << "Constraint not satisfied at step " << steps << ": K_eq*B - C = " << constraint_residual;

    time += dt;
    steps++;
  }

  // Verify constraint is satisfied: C = K_eq * B
  double expected_c = state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx];
  double final_residual =
      state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];

  EXPECT_NEAR(state.variables_[0][C_idx], expected_c, 1.0e-6);
  EXPECT_NEAR(final_residual, 0.0, 1.0e-6);

  // Verify mass conservation: A + B should be conserved (approximately)
  double total = state.variables_[0][A_idx] + state.variables_[0][B_idx];
  EXPECT_NEAR(total, 1.0, 0.01);
}

/// @brief Test DAE solve with constraints and state reordering enabled
TEST(EquilibriumIntegration, DAESolveWithConstraintAndReorderState)
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  double k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  double k_eq = 2.0;
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(true)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");
  std::size_t b_c_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 400.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly
  solver.UpdateStateParameters(state);
  double expected_k_eq = ComputeEquilibriumConstant(k_eq, -2400.0, 400.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")], expected_k_eq, 1e-10);

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Reordered DAE solve did not converge at time=" << time;

    // Constraint should hold at each step
    double residual = state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
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
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");
  auto d = Species("D");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C, D } };

  double k = 1.0;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  double k_eq1 = 3.0;
  double k_eq2 = 5.0;
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      C,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { C, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq1, .delta_H = -2400.0 }));
  constraints.push_back(EquilibriumConstraint(
      "B_D_eq",
      D,
      std::vector<StoichSpecies>{ { B, 1.0 } },
      std::vector<StoichSpecies>{ { D, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq2, .delta_H = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");
  std::size_t d_idx = state.variable_map_.at("D");
  std::size_t b_c_eq_idx = state.custom_rate_parameter_map_.at("B_C_eq");
  std::size_t b_d_eq_idx = state.custom_rate_parameter_map_.at("B_D_eq");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.variables_[0][D_idx] = 0.0;
  state.conditions_[0].temperature_ = 300.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly for both constraints
  solver.UpdateStateParameters(state);
  double expected_k_eq1 = ComputeEquilibriumConstant(k_eq1, -2400.0, 300.0);
  double expected_k_eq2 = ComputeEquilibriumConstant(k_eq2, -2400.0, 300.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_C_eq_idx], expected_k_eq1, 1e-10);
  EXPECT_NEAR(state.custom_rate_parameters_[0][B_D_eq_idx], expected_k_eq2, 1e-10);

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Coupled constraints did not converge at time=" << time;

    double residual1 =
        state.custom_rate_parameters_[0][B_C_eq_idx] * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    double residual2 =
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
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Reaction: C -> A (to produce A for the equilibrium)
  double k = 0.5;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ C })
                    .SetProducts({ { A, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Equilibrium constraint: K_eq * [A]^2 - [B] = 0
  double k_eq = 10.0;
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A2_B_eq",
      B,
      std::vector<StoichSpecies>{ { A, 2.0 } },
      std::vector<StoichSpecies>{ { B, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 }));

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");
  std::size_t a2_b_eq_idx = state.custom_rate_parameter_map_.at("A2_B_eq");

  // Start with some A so the constraint has something to work with
  state.variables_[0][A_idx] = 0.1;
  state.variables_[0][B_idx] = k_eq * 0.1 * 0.1;  // B = K_eq * A^2
  state.variables_[0][C_idx] = 1.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify UpdateStateParameters calculates K_eq correctly
  solver.UpdateStateParameters(state);
  double expected_k_eq = ComputeEquilibriumConstant(k_eq, -2400.0, 298.0);
  EXPECT_NEAR(state.custom_rate_parameters_[0][A2_B_eq_idx], expected_k_eq, 1e-10);

  double dt = 0.001;
  double total_time = 0.05;
  double time = 0.0;

  while (time < total_time)
  {
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "NonUnit stoich did not converge at time=" << time;

    // Constraint: K_eq * [A]^2 - [B] = 0
    double a_val = state.variables_[0][A_idx];
    double b_val = state.variables_[0][B_idx];
    double residual = state.custom_rate_parameters_[0][A2_B_eq_idx] * a_val * a_val - b_val;
    EXPECT_NEAR(residual, 0.0, 1.0e-5) << "Constraint violated at time=" << time;

    time += dt;
  }
}