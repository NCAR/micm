// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/equilibrium_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

/// @brief Test ConstraintSet API directly (unit-level test for DAE infrastructure)
///        Uses replace-state-rows mode: AB's row (index 2) is replaced by the constraint.
TEST(EquilibriumIntegration, ConstraintSetAPITest)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto AB = micm::Species("AB");

  // Create constraint: A + B <-> AB with K_eq = 1000
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "A_B_eq",
      std::vector<micm::StoichSpecies>{ { A, 1.0 }, { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { AB, 1.0 } },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 }
  };

  micm::ConstraintSet set(std::move(constraints), variable_map);

  // Test at equilibrium: [A] = 0.0312, [B] = 0.0312, [AB] = 0.9737
  // K_eq * [A] * [B] = 1000 * 0.0312^2 ≈ 0.9734
  // Should approximately equal [AB] = 0.9737
  micm::Matrix<double> state(1, 3);
  state[0][0] = 0.0312;  // A
  state[0][1] = 0.0312;  // B
  state[0][2] = 0.9737;  // AB - calculated to be at equilibrium

  micm::Matrix<double> forcing(1, 3, 0.0);
  set.AddForcingTerms(state, forcing);

  // At equilibrium, residual replaces AB's row (index 2)
  // G = K_eq * [A] * [B] - [AB] = 1000 * 0.0312 * 0.0312 - 0.9737 ≈ 0
  EXPECT_NEAR(forcing[0][2], 0.0, 0.01);
}

/// @brief Test SetConstraints API integration - verifies solver builds and runs with constraints
///
/// This test verifies that the SolverBuilder correctly accepts constraints via SetConstraints
/// and that the resulting solver can be built without errors. The actual DAE solving
/// requires additional solver infrastructure for proper algebraic equation handling.
TEST(EquilibriumIntegration, SetConstraintsAPIWorks)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B
  double k = 0.1;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // Create an equilibrium constraint
  // C is an algebraic variable (not in any kinetic reaction)
  double K_eq = 10.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  // Build solver with constraints - this verifies the API works
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify constraint metadata
  ASSERT_EQ(state.state_size_, 3);       // A, B, and C
  ASSERT_EQ(state.constraint_size_, 1);
  ASSERT_TRUE(state.variable_map_.count("A") > 0);
  ASSERT_TRUE(state.variable_map_.count("B") > 0);
  ASSERT_TRUE(state.variable_map_.count("C") > 0);

  // Verify mass-matrix diagonal.
  // Size is species only; constrained species rows are algebraic (0), others are ODE (1).
  ASSERT_EQ(state.upper_left_identity_diagonal_.size(), 3);
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("A")], 1.0);  // A is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("B")], 1.0);  // B is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("C")], 0.0);  // C is algebraic

  // Verify state can be initialized
  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.1;
  state.variables_[0][C_idx] = K_eq * 0.1;  // C = K_eq * B
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify CalculateRateConstants works
  solver.CalculateRateConstants(state);
}

/// @brief Test SetConstraints API with multiple constraints
///
/// Verifies that multiple constraints can be added via SetConstraints and the solver
/// builds correctly with the proper state structure.
TEST(EquilibriumIntegration, SetConstraintsAPIMultipleConstraints)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  auto D = micm::Species("D");
  auto E = micm::Species("E");
  auto F = micm::Species("F");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C, D, E, F } };

  // Simple kinetic reactions
  micm::Process rxn1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ A })
                           .SetProducts({ { B, 1 } })
                           .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.5, .B_ = 0, .C_ = 0 }))
                           .SetPhase(gas_phase)
                           .Build();

  micm::Process rxn2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ D })
                           .SetProducts({ { E, 1 } })
                           .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 0.2, .B_ = 0, .C_ = 0 }))
                           .SetPhase(gas_phase)
                           .Build();

  // Two equilibrium constraints
  double K_eq1 = 5.0;
  double K_eq2 = 20.0;

  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq1));
  constraints.push_back(micm::EquilibriumConstraint(
      "E_F_eq",
      std::vector<micm::StoichSpecies>{ { E, 1.0 } },
      std::vector<micm::StoichSpecies>{ { F, 1.0 } },
      K_eq2));

  // Build solver with multiple constraints
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn1, rxn2 })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify state structure
  ASSERT_EQ(state.state_size_, 6);       // A, B, C, D, E, F
  ASSERT_EQ(state.constraint_size_, 2);  // Two constraints
  ASSERT_TRUE(state.variable_map_.count("A") > 0);
  ASSERT_TRUE(state.variable_map_.count("B") > 0);
  ASSERT_TRUE(state.variable_map_.count("C") > 0);
  ASSERT_TRUE(state.variable_map_.count("D") > 0);
  ASSERT_TRUE(state.variable_map_.count("E") > 0);
  ASSERT_TRUE(state.variable_map_.count("F") > 0);

  // Verify mass-matrix diagonal.
  // Size is species only; constrained species rows are algebraic (0), others are ODE (1).
  ASSERT_EQ(state.upper_left_identity_diagonal_.size(), 6);
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("A")], 1.0);  // A is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("B")], 1.0);  // B is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("C")], 0.0);  // C is algebraic
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("D")], 1.0);  // D is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("E")], 1.0);  // E is ODE
  EXPECT_EQ(state.upper_left_identity_diagonal_[state.variable_map_.at("F")], 0.0);  // F is algebraic

  // Initialize state and verify CalculateRateConstants works
  state.variables_[0][state.variable_map_.at("A")] = 1.0;
  state.variables_[0][state.variable_map_.at("B")] = 0.1;
  state.variables_[0][state.variable_map_.at("C")] = K_eq1 * 0.1;  // C should satisfy C = K_eq1 * B
  state.variables_[0][state.variable_map_.at("D")] = 0.5;
  state.variables_[0][state.variable_map_.at("E")] = 0.05;
  state.variables_[0][state.variable_map_.at("F")] = K_eq2 * 0.05;  // F should satisfy F = K_eq2 * E
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  solver.CalculateRateConstants(state);
}

/// @brief Test DAE solving - actually calls Solve() with algebraic constraints
/// This exercises mass-matrix DAE enforcement where constrained species rows are algebraic.
TEST(EquilibriumIntegration, DAESolveWithConstraint)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  // This couples B (ODE variable) to C (algebraic variable)
  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Initial conditions: A=1, B=0, C=0
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Solve with smaller time steps
  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;
  int steps = 0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);

    if (result.state_ != micm::SolverState::Converged)
    {
      FAIL() << "DAE solve did not converge at step " << steps << ", time=" << time;
    }

    // Verify constraint is maintained by the solver
    double constraint_residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(constraint_residual, 0.0, 1.0e-6)
        << "Constraint not satisfied at step " << steps 
        << ": K_eq*B - C = " << constraint_residual;

    time += dt;
    steps++;
  }

  // Verify constraint is satisfied: C = K_eq * B
  double expected_C = K_eq * state.variables_[0][B_idx];
  double final_residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];

  EXPECT_NEAR(state.variables_[0][C_idx], expected_C, 1.0e-6);
  EXPECT_NEAR(final_residual, 0.0, 1.0e-6);

  // Verify mass conservation: A + B should be conserved (approximately)
  double total = state.variables_[0][A_idx] + state.variables_[0][B_idx];
  EXPECT_NEAR(total, 1.0, 0.01);
}

/// @brief Test DAE solve with constraints and state reordering enabled
TEST(EquilibriumIntegration, DAESolveWithConstraintAndReorderState)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  // Simple reaction: A -> B with rate k
  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(true)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Verify algebraic/ODE split is preserved under reordering
  EXPECT_EQ(state.upper_left_identity_diagonal_[A_idx], 1.0);
  EXPECT_EQ(state.upper_left_identity_diagonal_[B_idx], 1.0);
  EXPECT_EQ(state.upper_left_identity_diagonal_[C_idx], 0.0);

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged)
        << "Reordered DAE solve did not converge at time=" << time;

    // Constraint should hold at each step
    double residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-6) << "Constraint violated at time=" << time;

    time += dt;
  }

  // Mass conservation for A -> B
  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test DAE solve with FourStageDifferentialAlgebraicRosenbrockParameters
TEST(EquilibriumIntegration, DAESolveWithFourStageDAEParameters)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "FourStageDAE did not converge at time=" << time;

    double residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-6) << "Constraint violated at time=" << time;

    time += dt;
  }

  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test DAE solve with SixStageDifferentialAlgebraicRosenbrockParameters
TEST(EquilibriumIntegration, DAESolveWithSixStageDAEParameters)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "SixStageDAE did not converge at time=" << time;

    double residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-6) << "Constraint violated at time=" << time;

    time += dt;
  }

  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test DAE solve with two coupled constraints sharing a species
/// A -> B (kinetic), B <-> C (constraint 1), B <-> D (constraint 2)
TEST(EquilibriumIntegration, DAESolveWithTwoCoupledConstraints)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  auto D = micm::Species("D");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C, D } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq1 = 3.0;
  double K_eq2 = 5.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq1));
  constraints.push_back(micm::EquilibriumConstraint(
      "B_D_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { D, 1.0 } },
      K_eq2));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");
  std::size_t D_idx = state.variable_map_.at("D");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.variables_[0][D_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "Coupled constraints did not converge at time=" << time;

    double residual1 = K_eq1 * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    double residual2 = K_eq2 * state.variables_[0][B_idx] - state.variables_[0][D_idx];
    EXPECT_NEAR(residual1, 0.0, 1.0e-6) << "Constraint 1 violated at time=" << time;
    EXPECT_NEAR(residual2, 0.0, 1.0e-6) << "Constraint 2 violated at time=" << time;

    time += dt;
  }

  // Both constraints should be satisfied at the end
  EXPECT_NEAR(K_eq1 * state.variables_[0][B_idx], state.variables_[0][C_idx], 1.0e-6);
  EXPECT_NEAR(K_eq2 * state.variables_[0][B_idx], state.variables_[0][D_idx], 1.0e-6);
  // Mass conservation: A + B should be conserved
  EXPECT_NEAR(state.variables_[0][A_idx] + state.variables_[0][B_idx], 1.0, 0.01);
}

/// @brief Test conservation law: A + B + C/K_eq should be conserved
TEST(EquilibriumIntegration, DAEConservationLaw)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Start at equilibrium for the constraint
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.1;
  state.variables_[0][C_idx] = K_eq * 0.1;  // C = K_eq * B at equilibrium
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // The kinetic reaction A -> B doesn't create or destroy C directly.
  // Since C = K_eq * B, the "effective B" including its constrained partner is B + C/K_eq = 2*B.
  // So A + B is the conserved quantity for the kinetic system.
  double initial_A_plus_B = state.variables_[0][A_idx] + state.variables_[0][B_idx];

  double dt = 0.001;
  double total_time = 0.5;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged);
    time += dt;
  }

  double final_A_plus_B = state.variables_[0][A_idx] + state.variables_[0][B_idx];
  EXPECT_NEAR(final_A_plus_B, initial_A_plus_B, 0.01);

  // Constraint should still be satisfied
  EXPECT_NEAR(K_eq * state.variables_[0][B_idx], state.variables_[0][C_idx], 1.0e-6);
}

/// @brief Test that the solver handles large K_eq (stiff coupling) without NaN
TEST(EquilibriumIntegration, DAESolveStiffCoupling)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // Large equilibrium constant creates stiff coupling
  double K_eq = 1000.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.05;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);

    // Should not produce NaN or Inf
    ASSERT_FALSE(std::isnan(state.variables_[0][A_idx])) << "NaN in A at time=" << time;
    ASSERT_FALSE(std::isnan(state.variables_[0][B_idx])) << "NaN in B at time=" << time;
    ASSERT_FALSE(std::isnan(state.variables_[0][C_idx])) << "NaN in C at time=" << time;
    ASSERT_FALSE(std::isinf(state.variables_[0][A_idx])) << "Inf in A at time=" << time;
    ASSERT_FALSE(std::isinf(state.variables_[0][B_idx])) << "Inf in B at time=" << time;
    ASSERT_FALSE(std::isinf(state.variables_[0][C_idx])) << "Inf in C at time=" << time;

    // Solver should converge (or at least not crash)
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "Stiff DAE did not converge at time=" << time;

    double residual = K_eq * state.variables_[0][B_idx] - state.variables_[0][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-4) << "Constraint violated at time=" << time;

    time += dt;
  }
}

/// @brief Test DAE solve with non-unit stoichiometric coefficient
/// 2A <-> B means K_eq * [A]^2 - [B] = 0
TEST(EquilibriumIntegration, DAESolveWithNonUnitStoichiometry)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  // Reaction: C -> A (to produce A for the equilibrium)
  double k = 0.5;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ C })
                          .SetProducts({ { A, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // Equilibrium constraint: K_eq * [A]^2 - [B] = 0
  double K_eq = 10.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "A2_B_eq",
      std::vector<micm::StoichSpecies>{ { A, 2.0 } },
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Start with some A so the constraint has something to work with
  state.variables_[0][A_idx] = 0.1;
  state.variables_[0][B_idx] = K_eq * 0.1 * 0.1;  // B = K_eq * A^2
  state.variables_[0][C_idx] = 1.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  double dt = 0.001;
  double total_time = 0.05;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "NonUnit stoich did not converge at time=" << time;

    // Constraint: K_eq * [A]^2 - [B] = 0
    double A_val = state.variables_[0][A_idx];
    double B_val = state.variables_[0][B_idx];
    double residual = K_eq * A_val * A_val - B_val;
    EXPECT_NEAR(residual, 0.0, 1.0e-5) << "Constraint violated at time=" << time;

    time += dt;
  }
}

/// @brief Test DAE solve with multiple grid cells
/// Exercises the per-cell iteration in AddForcingTerms, SubtractJacobianTerms,
/// AlphaMinusJacobian, and NormalizedError with different initial conditions per cell.
TEST(EquilibriumIntegration, DAESolveMultiGridCell)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  const std::size_t num_cells = 3;
  auto state = solver.GetState(num_cells);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Different initial conditions per cell, all starting at constraint equilibrium (C = K_eq * B)
  // Cell 0: high A, small B
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.01;
  state.variables_[0][C_idx] = K_eq * 0.01;

  // Cell 1: moderate A and B
  state.variables_[1][A_idx] = 0.5;
  state.variables_[1][B_idx] = 0.2;
  state.variables_[1][C_idx] = K_eq * 0.2;

  // Cell 2: low A, high B
  state.variables_[2][A_idx] = 0.1;
  state.variables_[2][B_idx] = 0.8;
  state.variables_[2][C_idx] = K_eq * 0.8;

  for (std::size_t cell = 0; cell < num_cells; ++cell)
  {
    state.conditions_[cell].temperature_ = 298.0;
    state.conditions_[cell].pressure_ = 101325.0;
  }

  // Record initial A+B per cell for conservation check
  std::vector<double> initial_A_plus_B(num_cells);
  for (std::size_t cell = 0; cell < num_cells; ++cell)
  {
    initial_A_plus_B[cell] = state.variables_[cell][A_idx] + state.variables_[cell][B_idx];
  }

  double dt = 0.001;
  double total_time = 0.1;
  double time = 0.0;

  while (time < total_time)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(dt, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged) << "Multi-grid DAE did not converge at time=" << time;

    // Check constraint satisfaction per cell
    for (std::size_t cell = 0; cell < num_cells; ++cell)
    {
      double residual = K_eq * state.variables_[cell][B_idx] - state.variables_[cell][C_idx];
      EXPECT_NEAR(residual, 0.0, 1.0e-6)
          << "Constraint violated in cell " << cell << " at time=" << time;

      // No NaN or Inf in any cell
      ASSERT_FALSE(std::isnan(state.variables_[cell][A_idx])) << "NaN in A, cell " << cell;
      ASSERT_FALSE(std::isnan(state.variables_[cell][B_idx])) << "NaN in B, cell " << cell;
      ASSERT_FALSE(std::isnan(state.variables_[cell][C_idx])) << "NaN in C, cell " << cell;
    }

    time += dt;
  }

  // Final checks per cell
  for (std::size_t cell = 0; cell < num_cells; ++cell)
  {
    // Constraint satisfied
    EXPECT_NEAR(K_eq * state.variables_[cell][B_idx], state.variables_[cell][C_idx], 1.0e-6)
        << "Final constraint violated in cell " << cell;

    // Mass conservation: A + B should be conserved
    double final_A_plus_B = state.variables_[cell][A_idx] + state.variables_[cell][B_idx];
    EXPECT_NEAR(final_A_plus_B, initial_A_plus_B[cell], 0.01)
        << "Mass not conserved in cell " << cell;

    // ODE variables non-negative
    EXPECT_GE(state.variables_[cell][A_idx], 0.0) << "Negative A in cell " << cell;
    EXPECT_GE(state.variables_[cell][B_idx], 0.0) << "Negative B in cell " << cell;
  }

  // Verify cells evolved differently (different initial conditions should give different results)
  EXPECT_NE(state.variables_[0][B_idx], state.variables_[1][B_idx]);
  EXPECT_NE(state.variables_[1][B_idx], state.variables_[2][B_idx]);
}

/// @brief Test that Max(0.0) clamping doesn't break algebraic variables
/// Algebraic variables should not be clamped since they're determined by constraints
TEST(EquilibriumIntegration, DAEClampingDoesNotBreakAlgebraicVariables)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Verify that constraints_replace_state_rows_ is set
  ASSERT_TRUE(state.constraints_replace_state_rows_);

  // After a solve, C should be K_eq * B, and the constraint should be satisfied
  // This verifies that clamping doesn't interfere with algebraic variable values
  solver.CalculateRateConstants(state);
  auto result = solver.Solve(0.01, state);
  ASSERT_EQ(result.state_, micm::SolverState::Converged);

  // B should be positive after one step (A is decaying into B)
  EXPECT_GT(state.variables_[0][B_idx], 0.0);
  // C should be K_eq * B
  EXPECT_NEAR(K_eq * state.variables_[0][B_idx], state.variables_[0][C_idx], 1.0e-6);
  // ODE variables (A, B) should be non-negative (clamped)
  EXPECT_GE(state.variables_[0][A_idx], 0.0);
  EXPECT_GE(state.variables_[0][B_idx], 0.0);
}

/// @brief Regression test for State copy preserving solver-specific temporary variables (Clone fix)
/// Copies a solver-generated State and solves with both original and copy to verify
/// that polymorphic Clone() correctly preserves RosenbrockTemporaryVariables.
TEST(EquilibriumIntegration, DAEStateCopyAndSolve)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.1;
  state.variables_[0][C_idx] = K_eq * 0.1;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Copy the state (exercises Clone())
  auto state_copy = state;

  // Solve with the original state
  solver.CalculateRateConstants(state);
  auto result1 = solver.Solve(0.01, state);
  ASSERT_EQ(result1.state_, micm::SolverState::Converged);

  // Solve with the copied state
  solver.CalculateRateConstants(state_copy);
  auto result2 = solver.Solve(0.01, state_copy);
  ASSERT_EQ(result2.state_, micm::SolverState::Converged);

  // Both should produce identical results
  EXPECT_DOUBLE_EQ(state.variables_[0][A_idx], state_copy.variables_[0][A_idx]);
  EXPECT_DOUBLE_EQ(state.variables_[0][B_idx], state_copy.variables_[0][B_idx]);
  EXPECT_DOUBLE_EQ(state.variables_[0][C_idx], state_copy.variables_[0][C_idx]);

  // Both should satisfy the constraint
  EXPECT_NEAR(K_eq * state.variables_[0][B_idx], state.variables_[0][C_idx], 1.0e-6);
  EXPECT_NEAR(K_eq * state_copy.variables_[0][B_idx], state_copy.variables_[0][C_idx], 1.0e-6);
}

/// @brief Regression test for duplicate dependency accumulation in Jacobian
/// When a species appears on both sides of a constraint (e.g., A + B <-> A + C),
/// the Jacobian entries for that species should accumulate correctly via -=.
TEST(EquilibriumIntegration, DAEOverlappingSpeciesJacobian)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  // Create constraint where A appears on both sides: A + B <-> A + C
  // G = K_eq * [A] * [B] - [A] * [C] = 0
  // At equilibrium: K_eq * [B] = [C] (when [A] > 0)
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "AB_AC_eq",
      std::vector<micm::StoichSpecies>{ { A, 1.0 }, { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { A, 1.0 }, { C, 1.0 } },
      2.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "C", 2 }
  };

  // Algebraic species is A (first product), so constraint replaces row 0
  micm::ConstraintSet set(std::move(constraints), variable_map);

  // Build a sparse Jacobian and set flat IDs
  auto nonzero = set.NonZeroJacobianElements();
  auto jacobian = micm::BuildJacobian<micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>>(
      nonzero, 1, 3, false);
  set.SetJacobianFlatIds(jacobian);

  // Set concentrations: A=1.0, B=0.5, C=1.0 (K_eq*B = 2*0.5 = 1.0 = C, at equilibrium)
  micm::Matrix<double> state(1, 3);
  state[0][0] = 1.0;   // A
  state[0][1] = 0.5;   // B
  state[0][2] = 1.0;   // C

  // Zero the Jacobian and subtract terms
  jacobian.Fill(0.0);
  set.SubtractJacobianTerms(state, jacobian);

  // Verify no NaN/Inf in the Jacobian
  for (std::size_t row = 0; row < 3; ++row)
  {
    for (std::size_t col = 0; col < 3; ++col)
    {
      if (jacobian.IsZero(row, col))
        continue;
      EXPECT_FALSE(std::isnan(jacobian[0][row][col]))
          << "NaN in Jacobian[" << row << "][" << col << "]";
      EXPECT_FALSE(std::isinf(jacobian[0][row][col]))
          << "Inf in Jacobian[" << row << "][" << col << "]";
    }
  }

  // Residual at equilibrium should be near zero
  micm::Matrix<double> forcing(1, 3, 0.0);
  set.AddForcingTerms(state, forcing);
  EXPECT_NEAR(forcing[0][0], 0.0, 1.0e-12);
}

/// @brief Regression test for constraint-only species with strict unused-species check
/// Species used only by constraints (not by any kinetic reaction) should not trigger
/// the unused-species error when SetIgnoreUnusedSpecies(false) is used.
TEST(EquilibriumIntegration, DAEConstraintOnlySpeciesNotUnused)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ A, B, C } };

  // Only reaction involves A and B; C is constraint-only
  double k = 1.0;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  // C is the algebraic species - only referenced by the constraint
  double K_eq = 2.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  // Build with SetIgnoreUnusedSpecies(false) - should NOT throw for constraint-only species
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  EXPECT_NO_THROW({
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetIgnoreUnusedSpecies(false)
                      .SetReorderState(false)
                      .Build();
  });
}
