// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "stub_aerosol_with_constraints.hpp"

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>

/// @brief Verify that AddExternalModel wraps both processes and constraints for a constrained model
TEST(ExternalModelConstraints, AddExternalModelWithConstraints)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify constraint metadata
  EXPECT_EQ(state.state_size_, 2);       // A_GAS, AEROSOL.A_AQ
  EXPECT_EQ(state.constraint_size_, 1);  // one algebraic constraint

  // Verify mass matrix diagonal: A_AQ row should be algebraic (0.0)
  auto i_gas = state.variable_map_.at("A_GAS");
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_gas], 1.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);
}

/// @brief Verify that AddExternalModel works for a process-only model (no constraints)
TEST(ExternalModelConstraints, AddExternalModelProcessOnly)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  // No total_mass → constraints disabled
  StubAerosolWithConstraints aerosol(0.01);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // No constraints — pure ODE mode
  EXPECT_EQ(state.constraint_size_, 0);
  auto i_gas = state.variable_map_.at("A_GAS");
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_gas], 1.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 1.0);
}

/// @brief Verify that the DAE solver enforces the mass conservation constraint
TEST(ExternalModelConstraints, DAESolveEnforcesConservation)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  double k = 0.1;
  StubAerosolWithConstraints aerosol(k, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Initialize: most mass in gas phase
  state.variables_[0][state.variable_map_.at("A_GAS")] = 0.9;
  state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")] = 0.1;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Solve several time steps and verify conservation
  double dt = 10.0;
  for (int step = 0; step < 20; ++step)
  {
    auto result = solver.Solve(dt, state);
    EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Step " << step;

    double sum = state.variables_[0][state.variable_map_.at("A_GAS")]
               + state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")];
    EXPECT_NEAR(sum, total, 1e-4) << "Conservation violated at step " << step;
  }
}

/// @brief Verify that external model constraints combine with built-in SetConstraints
TEST(ExternalModelConstraints, CombinedBuiltInAndExternalConstraints)
{
  auto A_GAS = micm::Species("A_GAS");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  micm::Phase gas_phase{ "gas", { A_GAS, B, C } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  // Built-in constraint: B <-> C equilibrium
  double K_eq = 5.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      K_eq));

  // Process: A_GAS -> B
  double k_rxn = 0.05;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A_GAS })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_rxn, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // 2 algebraic variables: C (from built-in) and AEROSOL.A_AQ (from external)
  EXPECT_EQ(state.constraint_size_, 2);

  // Verify mass matrix diagonal
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  auto i_c = state.variable_map_.at("C");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_c], 0.0);
}

/// @brief Verify backward-compatible AddExternalModelProcesses() still works
TEST(ExternalModelConstraints, BackwardCompatAddExternalModelProcesses)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  StubAerosolWithConstraints aerosol(0.01);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Use the old API — should still work (processes only, no constraints)
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModelProcesses(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  EXPECT_EQ(state.constraint_size_, 0);
}

/// @brief Verify AddExternalModelConstraints() adds only constraints (processes added separately)
TEST(ExternalModelConstraints, AddExternalModelConstraintsOnly)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModelProcesses(aerosol)    // processes only
                    .AddExternalModelConstraints(aerosol)   // constraints only
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Constraint should be active
  EXPECT_EQ(state.constraint_size_, 1);
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);

  // Solve and verify conservation
  state.variables_[0][state.variable_map_.at("A_GAS")] = 0.8;
  state.variables_[0][i_aq] = 0.2;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  auto result = solver.Solve(50.0, state);
  EXPECT_EQ(result.state_, micm::SolverState::Converged);

  double sum = state.variables_[0][state.variable_map_.at("A_GAS")]
             + state.variables_[0][i_aq];
  EXPECT_NEAR(sum, total, 1e-4);
}

/// @brief Verify multiple grid cells work with external model constraints
TEST(ExternalModelConstraints, MultiGridCell)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.1, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  const int num_cells = 3;
  auto state = solver.GetState(num_cells);

  // Different initial conditions per cell
  for (int c = 0; c < num_cells; ++c)
  {
    double gas_frac = 0.9 - 0.2 * c;
    state.variables_[c][state.variable_map_.at("A_GAS")] = gas_frac;
    state.variables_[c][state.variable_map_.at("AEROSOL.A_AQ")] = total - gas_frac;
    state.conditions_[c].temperature_ = 298.0;
    state.conditions_[c].pressure_ = 101325.0;
  }

  auto result = solver.Solve(50.0, state);
  EXPECT_EQ(result.state_, micm::SolverState::Converged);

  for (int c = 0; c < num_cells; ++c)
  {
    double sum = state.variables_[c][state.variable_map_.at("A_GAS")]
               + state.variables_[c][state.variable_map_.at("AEROSOL.A_AQ")];
    EXPECT_NEAR(sum, total, 1e-4) << "Conservation violated in cell " << c;
  }
}
