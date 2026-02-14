// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "stub_aerosol_1.hpp"
#include "stub_aerosol_2.hpp"

#include <micm/CPU.hpp>

#include <gtest/gtest.h>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <vector>

/// @brief  Helper function to create a system with two stub aerosol models
/// @return A tuple containing the system, the first stub aerosol model, the second stub aerosol model, and a map of phase names to Phase objects
inline std::tuple<micm::System, StubAerosolModel, AnotherStubAerosolModel, std::map<std::string, micm::Phase>> CreateSystemWithStubAerosolModels()
{
  // Create a simple chemical system
  auto foo = micm::Species("FO2"); // species that can partition to condensed phase
  auto bar = micm::Species("BAR"); // gas-phase only species
  auto baz = micm::Species("BAZ"); // condensed-phase only species
  auto qux = micm::Species("QUX"); // condensed-phase only species

  auto gas   = micm::Phase("GAS", std::vector<micm::PhaseSpecies>({ foo, bar }));        // gas phase
  auto quux  = micm::Phase("QUUX", std::vector<micm::PhaseSpecies>({ baz, qux }));       // condensed aerosol or cloud phase
  auto corge = micm::Phase("CORGE", std::vector<micm::PhaseSpecies>({ foo, baz, qux })); // another condensed aerosol or cloud phase
  std::map<std::string, micm::Phase> phases = {
    { "GAS", gas },
    { "QUUX", quux },
    { "CORGE", corge }
  };

  // Create instances of each stub aerosol model
  StubAerosolModel::RateConstants rate_constants = {
    .fo2_gas_to_mode2_corge = STUB1_RATE_CONSTANT_FO2_CORGE,
    .baz_mode1_to_mode2_quux = STUB1_RATE_CONSTANT_BAZ_QUUX
  };
  auto aerosol_1 = StubAerosolModel("STUB1", std::vector<micm::Phase>({ quux, corge }), rate_constants);
  auto aerosol_2 = AnotherStubAerosolModel("STUB2", std::vector<micm::Phase>({ quux, corge }));

  // Create a system containing the gas phase and both aerosol models
  auto system = micm::System({
    .gas_phase_ = gas,
    .external_models_ = { aerosol_1, aerosol_2 }
  });

  return { system, aerosol_1, aerosol_2, phases };
}

/// @brief Test that state includes stub aerosol model variables
template<class BuilderPolicy>
void test_state_includes_stub_aerosol_model(BuilderPolicy builder)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto solver = builder.SetSystem(system)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and ensure that the size and labels match expectations
  auto state = solver.GetState();
  EXPECT_EQ(
      state.variable_map_.size(),
      system.gas_phase_.UniqueNames().size() + aerosol_1.StateVariableNames().size() + aerosol_2.StateVariableNames().size());
  EXPECT_EQ(std::get<0>(aerosol_1.StateSize()), aerosol_1.StateVariableNames().size());
  EXPECT_EQ(std::get<0>(aerosol_2.StateSize()), aerosol_2.StateVariableNames().size());
  EXPECT_EQ(std::get<1>(aerosol_1.StateSize()), aerosol_1.StateParameterNames().size());
  EXPECT_EQ(std::get<1>(aerosol_2.StateSize()), aerosol_2.StateParameterNames().size());

  // Assemble the full list of expected variable names
  std::vector<std::string> expected_names;
  auto gas_names = system.gas_phase_.SpeciesNames();
  expected_names.insert(expected_names.end(), gas_names.begin(), gas_names.end());
  auto aerosol1_names = aerosol_1.StateVariableNames();
  expected_names.insert(expected_names.end(), aerosol1_names.begin(), aerosol1_names.end());
  auto aerosol2_names = aerosol_2.StateVariableNames();
  expected_names.insert(expected_names.end(), aerosol2_names.begin(), aerosol2_names.end());

  // Ensure that each unique name exists in the variable map
  for (const auto& name : expected_names)
  {
    EXPECT_TRUE(state.variable_map_.find(name) != state.variable_map_.end());
  }

  // Check a few specific variable names to ensure correct naming convention
  EXPECT_TRUE(state.variable_map_.find("FO2") != state.variable_map_.end());
  EXPECT_TRUE(state.variable_map_.find("STUB1.MODE1.QUUX.QUX") != state.variable_map_.end());
  EXPECT_TRUE(state.variable_map_.find("STUB1.MODE2.CORGE.FO2") != state.variable_map_.end());
  EXPECT_TRUE(state.variable_map_.find("STUB2.MODE1.NUMBER") != state.variable_map_.end());
  EXPECT_TRUE(state.variable_map_.find("STUB2.MODE2.CORGE.BAZ") != state.variable_map_.end());
}

/// @brief Test updating state with stub aerosol model (single grid cell)
template<class BuilderPolicy>
void test_update_state_with_stub_aerosol_model(BuilderPolicy builder)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto solver = builder.SetSystem(system)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and set some values
  auto state = solver.GetState();
 
  // Set some gas-phase species by name
  state["FO2"] = 1.34;

  // Set some gas-phase species by species object
  auto bar = micm::Species("BAR");
  state[bar] = 2.53;

  // Set some condensed-phase species in the first aerosol model by unique name
  state["STUB1.MODE1.QUUX.BAZ"] = 0.75;
  state["STUB1.MODE2.CORGE.FO2"] = 1.23;

  // Set some condensed-phase species in the second aerosol model by unique name and species object
  auto& corge = phases.at("CORGE");
  auto qux = micm::Species("QUX");
  state[aerosol_2.Species(2, corge, qux)] = 0.42;
  // Set some condensed-phase species in the second aerosol model by index
  auto corge_baz_it = state.variable_map_.find("STUB2.MODE3.CORGE.BAZ");
  ASSERT_NE(corge_baz_it, state.variable_map_.end());
  auto corge_baz_index = corge_baz_it->second;
  state[corge_baz_index] = 0.33;

  // Set some aerosol-model-specific variables by unique name
  state["STUB2.MODE1.NUMBER"] = 1000.0;

  // Set some aerosol-model-specific variables using the model instance
  state[aerosol_2.Number(1)] = 500.0;

  // Verify that all values were set correctly
  EXPECT_DOUBLE_EQ(state["FO2"], 1.34);
  EXPECT_DOUBLE_EQ(state["BAR"], 2.53);
  EXPECT_DOUBLE_EQ(state["STUB1.MODE1.QUUX.BAZ"], 0.75);
  EXPECT_DOUBLE_EQ(state["STUB1.MODE2.CORGE.FO2"], 1.23);
  EXPECT_DOUBLE_EQ(state["STUB2.MODE3.CORGE.QUX"], 0.42);
  EXPECT_DOUBLE_EQ(state["STUB2.MODE3.CORGE.BAZ"], 0.33);
  EXPECT_DOUBLE_EQ(state["STUB2.MODE1.NUMBER"], 1000.0);
  EXPECT_DOUBLE_EQ(state["STUB2.MODE2.NUMBER"], 500.0);
  
  // Re-verify using indices from the variable map
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("FO2")->second], 1.34);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("BAR")->second], 2.53);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB1.MODE1.QUUX.BAZ")->second], 0.75);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB1.MODE2.CORGE.FO2")->second], 1.23);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB2.MODE3.CORGE.QUX")->second], 0.42);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB2.MODE3.CORGE.BAZ")->second], 0.33);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB2.MODE1.NUMBER")->second], 1000.0);
  EXPECT_DOUBLE_EQ(state[state.variable_map_.find("STUB2.MODE2.NUMBER")->second], 500.0);

  // Re-re-verify using species objects where applicable
  EXPECT_DOUBLE_EQ(state[micm::Species("FO2")], 1.34);
  EXPECT_DOUBLE_EQ(state[bar], 2.53);
  EXPECT_DOUBLE_EQ(state[aerosol_1.Species(0, phases["QUUX"], micm::Species("BAZ"))], 0.75);
  EXPECT_DOUBLE_EQ(state[aerosol_1.Species(1, phases["CORGE"], micm::Species("FO2"))], 1.23);
  EXPECT_DOUBLE_EQ(state[aerosol_2.Species(2, phases["CORGE"], micm::Species("QUX"))], 0.42);
  EXPECT_DOUBLE_EQ(state[aerosol_2.Species(2, phases["CORGE"], micm::Species("BAZ"))], 0.33);
  EXPECT_DOUBLE_EQ(state[aerosol_2.Number(0)], 1000.0);
  EXPECT_DOUBLE_EQ(state[aerosol_2.Number(1)], 500.0);
}

/// @brief Test updating multi-cell state with stub aerosol model
template<class BuilderPolicy>
void test_update_multi_cell_state_with_stub_aerosol_model(BuilderPolicy builder)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto solver = builder.SetSystem(system)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  const std::size_t num_cells = 3;

  // Get a state and set some values
  auto state = solver.GetState(num_cells);
 
  // Set some gas-phase species by name
  state["FO2"] = std::vector{ 1.34, 1.35, 1.36 };

  // Set some gas-phase species by species object
  auto bar = micm::Species("BAR");
  state[bar] = std::vector{ 2.53, 2.54, 2.55 };

  // Set some condensed-phase species in the first aerosol model by unique name
  state["STUB1.MODE1.QUUX.BAZ"] = std::vector{ 0.75, 0.76, 0.77 };
  state["STUB1.MODE2.CORGE.FO2"] = std::vector{ 1.23, 1.24, 1.25 };

  // Set some condensed-phase species in the second aerosol model by unique name and species object
  auto corge_it = phases.find("CORGE");
  ASSERT_NE(corge_it, phases.end());
  auto& corge = corge_it->second;
  auto qux = micm::Species("QUX");
  state[aerosol_2.Species(2, corge, qux)] = std::vector{ 0.42, 0.43, 0.44 };

  // Set some condensed-phase species in the second aerosol model by index
  auto corge_baz_it = state.variable_map_.find("STUB2.MODE3.CORGE.BAZ");
  ASSERT_NE(corge_baz_it, state.variable_map_.end());
  auto corge_baz_index = corge_baz_it->second;
  state[corge_baz_index] = std::vector{ 0.33, 0.34, 0.35 };

  // Set some aerosol-model-specific variables by unique name
  state["STUB2.MODE1.NUMBER"] = std::vector{ 1000.0, 1001.0, 1002.0 };

  // Set some aerosol-model-specific variables using the model instance
  state[aerosol_2.Number(1)] = std::vector{ 500.0, 501.0, 502.0 };

  // Verify that all values were set correctly
  EXPECT_EQ(state["FO2"], (std::vector{ 1.34, 1.35, 1.36 }));
  EXPECT_EQ(state["BAR"], (std::vector{ 2.53, 2.54, 2.55 }));
  EXPECT_EQ(state["STUB1.MODE1.QUUX.BAZ"], (std::vector{ 0.75, 0.76, 0.77 }));
  EXPECT_EQ(state["STUB1.MODE2.CORGE.FO2"], (std::vector{ 1.23, 1.24, 1.25 }));
  EXPECT_EQ(state["STUB2.MODE3.CORGE.QUX"], (std::vector{ 0.42, 0.43, 0.44 }));
  EXPECT_EQ(state["STUB2.MODE3.CORGE.BAZ"], (std::vector{ 0.33, 0.34, 0.35 }));
  EXPECT_EQ(state["STUB2.MODE1.NUMBER"], (std::vector{ 1000.0, 1001.0, 1002.0 }));
  EXPECT_EQ(state["STUB2.MODE2.NUMBER"], (std::vector{ 500.0, 501.0, 502.0 }));
  
  // Re-verify using indices from the variable map
  EXPECT_EQ(state[state.variable_map_.find("FO2")->second], (std::vector{ 1.34, 1.35, 1.36 }));
  EXPECT_EQ(state[state.variable_map_.find("BAR")->second], (std::vector{ 2.53, 2.54, 2.55 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB1.MODE1.QUUX.BAZ")->second], (std::vector{ 0.75, 0.76, 0.77 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB1.MODE2.CORGE.FO2")->second], (std::vector{ 1.23, 1.24, 1.25 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB2.MODE3.CORGE.QUX")->second], (std::vector{ 0.42, 0.43, 0.44 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB2.MODE3.CORGE.BAZ")->second], (std::vector{ 0.33, 0.34, 0.35 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB2.MODE1.NUMBER")->second], (std::vector{ 1000.0, 1001.0, 1002.0 }));
  EXPECT_EQ(state[state.variable_map_.find("STUB2.MODE2.NUMBER")->second], (std::vector{ 500.0, 501.0, 502.0 }));

  // Re-re-verify using species objects where applicable
  EXPECT_EQ(state[micm::Species("FO2")], (std::vector{ 1.34, 1.35, 1.36 }));
  EXPECT_EQ(state[bar], (std::vector{ 2.53, 2.54, 2.55 }));
  EXPECT_EQ(state[aerosol_1.Species(0, phases["QUUX"], micm::Species("BAZ"))], (std::vector{ 0.75, 0.76, 0.77 }));
  EXPECT_EQ(state[aerosol_1.Species(1, phases["CORGE"], micm::Species("FO2"))], (std::vector{ 1.23, 1.24, 1.25 }));
  EXPECT_EQ(state[aerosol_2.Species(2, phases["CORGE"], micm::Species("QUX"))], (std::vector{ 0.42, 0.43, 0.44 }));
  EXPECT_EQ(state[aerosol_2.Species(2, phases["CORGE"], micm::Species("BAZ"))], (std::vector{ 0.33, 0.34, 0.35 }));
  EXPECT_EQ(state[aerosol_2.Number(0)], (std::vector{ 1000.0, 1001.0, 1002.0 }));
  EXPECT_EQ(state[aerosol_2.Number(1)], (std::vector{ 500.0, 501.0, 502.0 }));
}

/// @brief Test single grid cell forcing calculation with stub aerosol model
template<class BuilderPolicy>
void test_single_cell_forcing_with_stub_aerosol_model(BuilderPolicy builder)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1)
                       .AddExternalModelProcesses(aerosol_2)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and set some values
  auto state = solver.GetState();
  state["FO2"] = 1.0;
  state["BAR"] = 2.0;
  state["STUB1.MODE1.QUUX.BAZ"] = 0.5;
  state["STUB1.MODE2.CORGE.FO2"] = 0.8;
  state["STUB2.MODE3.CORGE.QUX"] = 0.3;
  state["STUB2.MODE3.CORGE.BAZ"] = 0.2;
  state["STUB2.MODE1.NUMBER"] = 1000.0;
  state["STUB2.MODE2.NUMBER"] = 500.0;

  // Calculate forcing terms using the first aerosol model's forcing function
  using DenseMatrixPolicyType = decltype(state.variables_);
  auto forcing_function_1 = aerosol_1.ForcingFunction<DenseMatrixPolicyType>(state.custom_rate_parameter_map_, state.variable_map_);
  auto forcing_1 = state.variables_; // make a forcing matrix of the same size as the state variable matrix
  forcing_1 = 0.0; // initialize forcing terms to zero before calculation
  forcing_function_1(state.custom_rate_parameters_, state.variables_, forcing_1);

  // For the FO2 gas to mode 2 CORGE partitioning, we expect a loss of FO2 in the gas phase and a corresponding gain of FO2 in the mode 2 CORGE phase,
  // with values equal to the rate constant multiplied by the FO2 concentration
  auto fo2_gas_index_it = state.variable_map_.find("FO2");
  auto fo2_mode2_index_it = state.variable_map_.find("STUB1.MODE2.CORGE.FO2");
  ASSERT_NE(fo2_gas_index_it, state.variable_map_.end());
  ASSERT_NE(fo2_mode2_index_it, state.variable_map_.end());
  std::size_t fo2_gas_index = fo2_gas_index_it->second;
  std::size_t fo2_mode2_index = fo2_mode2_index_it->second;
  double expected_fo2_loss = -STUB1_RATE_CONSTANT_FO2_CORGE * state["FO2"];
  double expected_fo2_gain = -expected_fo2_loss; // should be equal and opposite to the loss
  EXPECT_DOUBLE_EQ(forcing_1[0][fo2_gas_index], expected_fo2_loss);
  EXPECT_DOUBLE_EQ(forcing_1[0][fo2_mode2_index], expected_fo2_gain);

  // For the baz mode 1 to mode 2 QUUX conversion, we expect a loss of baz in mode 1 QUUX and a corresponding gain of baz in mode 2 QUUX,
  // with values equal to the rate constant multiplied by the baz concentration in mode 1 QUUX
  auto baz_mode1_index_it = state.variable_map_.find("STUB1.MODE1.QUUX.BAZ");
  auto baz_mode2_index_it = state.variable_map_.find("STUB1.MODE2.QUUX.BAZ");
  ASSERT_NE(baz_mode1_index_it, state.variable_map_.end());
  ASSERT_NE(baz_mode2_index_it, state.variable_map_.end());
  std::size_t baz_mode1_index = baz_mode1_index_it->second;
  std::size_t baz_mode2_index = baz_mode2_index_it->second;
  double expected_baz_loss = -STUB1_RATE_CONSTANT_BAZ_QUUX * state["STUB1.MODE1.QUUX.BAZ"];
  double expected_baz_gain = -expected_baz_loss; // should be equal and opposite to the loss
  EXPECT_DOUBLE_EQ(forcing_1[0][baz_mode1_index], expected_baz_loss);
  EXPECT_DOUBLE_EQ(forcing_1[0][baz_mode2_index], expected_baz_gain);
}

/// @brief Test single grid cell Jacobian calculation with stub aerosol model
template<class BuilderPolicy>
void test_single_cell_jacobian_with_stub_aerosol_model(BuilderPolicy builder)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1)
                       .AddExternalModelProcesses(aerosol_2)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and set some values
  auto state = solver.GetState();
  state["FO2"] = 1.0;
  state["BAR"] = 2.0;
  state["STUB1.MODE1.QUUX.BAZ"] = 0.5;
  state["STUB1.MODE2.CORGE.FO2"] = 0.8;
  state["STUB2.MODE3.CORGE.QUX"] = 0.3;
  state["STUB2.MODE3.CORGE.BAZ"] = 0.2;
  state["STUB2.MODE1.NUMBER"] = 1000.0;
  state["STUB2.MODE2.NUMBER"] = 500.0;

  // Calculate Jacobian terms using the first aerosol model's Jacobian function
  auto jacobian_1 = state.jacobian_; // make a Jacobian matrix of the same size as the system Jacobian
  jacobian_1 = 0.0; // initialize Jacobian terms to zero before calculation
  using DenseMatrixPolicyType = decltype(state.variables_);
  using SparseMatrixPolicyType = decltype(state.jacobian_);
  auto jacobian_function_1 = aerosol_1.JacobianFunction<DenseMatrixPolicyType, SparseMatrixPolicyType>(state.custom_rate_parameter_map_, state.variable_map_, jacobian_1);
  jacobian_function_1(state.custom_rate_parameters_, state.variables_, jacobian_1);

  // For the FO2 gas to mode 2 CORGE partitioning, we expect two non-zero Jacobian elements: a negative value equal to the rate constant in the column corresponding to FO2 and row corresponding to FO2 (representing the partial derivative of the FO2 loss with respect to FO2), and a positive value equal to the rate constant in the column corresponding to FO2 and row corresponding to mode 2 CORGE FO2 (representing the partial derivative of the FO2 gain with respect to FO2)
  auto fo2_gas_index_it = state.variable_map_.find("FO2");
  auto fo2_mode2_index_it = state.variable_map_.find("STUB1.MODE2.CORGE.FO2");
  ASSERT_NE(fo2_gas_index_it, state.variable_map_.end());
  ASSERT_NE(fo2_mode2_index_it, state.variable_map_.end());
  std::size_t fo2_gas_index = fo2_gas_index_it->second;
  std::size_t fo2_mode2_index = fo2_mode2_index_it->second;
  double expected_fo2_gas_partial = -STUB1_RATE_CONSTANT_FO2_CORGE;
  double expected_fo2_mode2_partial = STUB1_RATE_CONSTANT_FO2_CORGE;
  // Jacobian function calculates the negative of the partial derivatives for use in implicit solvers
  EXPECT_DOUBLE_EQ(jacobian_1[0][fo2_gas_index][fo2_gas_index], -expected_fo2_gas_partial);
  EXPECT_DOUBLE_EQ(jacobian_1[0][fo2_mode2_index][fo2_gas_index], -expected_fo2_mode2_partial);

  // For the baz mode 1 to mode 2 QUUX conversion, we expect two non-zero Jacobian elements: a negative value equal to the rate constant in the column corresponding to baz mode 1 QUUX and row corresponding to baz mode 1 QUUX (representing the partial derivative of the baz loss with respect to baz), and a positive value equal to the rate constant in the column corresponding to baz mode 1 QUUX and row corresponding to baz mode 2 QUUX (representing the partial derivative of the baz gain with respect to baz)
  auto baz_mode1_index_it = state.variable_map_.find("STUB1.MODE1.QUUX.BAZ");
  auto baz_mode2_index_it = state.variable_map_.find("STUB1.MODE2.QUUX.BAZ");
  ASSERT_NE(baz_mode1_index_it, state.variable_map_.end());
  ASSERT_NE(baz_mode2_index_it, state.variable_map_.end());
  std::size_t baz_mode1_index = baz_mode1_index_it->second;
  std::size_t baz_mode2_index = baz_mode2_index_it->second;
  double expected_baz_mode1_partial = -STUB1_RATE_CONSTANT_BAZ_QUUX;
  double expected_baz_mode2_partial = STUB1_RATE_CONSTANT_BAZ_QUUX;
  // Jacobian function calculates the negative of the partial derivatives for use in implicit solvers
  EXPECT_DOUBLE_EQ(jacobian_1[0][baz_mode1_index][baz_mode1_index], -expected_baz_mode1_partial);
  EXPECT_DOUBLE_EQ(jacobian_1[0][baz_mode2_index][baz_mode1_index], -expected_baz_mode2_partial);
}

/// @brief Test solving with stub aerosol models
template<class BuilderPolicy>
void test_solve_with_stub_aerosol_model_1(BuilderPolicy builder, double base_relative_tolerance = 5e-5)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1) // excluding aerosol 2 process for this test
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and set some initial values
  auto state = solver.GetState();

  double fo2_initial = 1.0;
  double baz_mode1_initial = 0.5;
  double fo2_mode2_initial = 0.8;
  state["FO2"] = fo2_initial;
  state["BAR"] = 2.0;
  state["STUB1.MODE1.QUUX.BAZ"] = baz_mode1_initial;
  state["STUB1.MODE2.CORGE.FO2"] = fo2_mode2_initial;
  state["STUB2.MODE3.CORGE.QUX"] = 0.3;
  state["STUB2.MODE3.CORGE.BAZ"] = 0.2;
  state["STUB2.MODE1.NUMBER"] = 1000.0;
  state["STUB2.MODE2.NUMBER"] = 500.0;
  
  // Calculate the analytical solution to verify the results
  double time_step = 10.0; // seconds
  double stub1_rxn1_delta = fo2_initial * (1.0 - std::exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
  double stub1_rxn2_delta = baz_mode1_initial * (1.0 - std::exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
  
  // Solve the system for a single time step
  solver.CalculateRateConstants(state);
  auto results = solver.Solve(time_step, state);

  // Make sure the solver reports success
  EXPECT_EQ(results.state_, micm::SolverState::Converged);

  // Verify that the state variables have been updated
  EXPECT_NEAR(state["FO2"], fo2_initial - stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  EXPECT_EQ(state["BAR"], 2.0);
  EXPECT_NEAR(state["STUB1.MODE1.QUUX.BAZ"], baz_mode1_initial - stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  EXPECT_NEAR(state["STUB1.MODE2.QUUX.BAZ"], stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  EXPECT_NEAR(state["STUB1.MODE2.CORGE.FO2"], fo2_mode2_initial + stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  EXPECT_EQ(state["STUB2.MODE3.CORGE.QUX"], 0.3);
  EXPECT_EQ(state["STUB2.MODE3.CORGE.BAZ"], 0.2);
  EXPECT_EQ(state["STUB2.MODE1.NUMBER"], 1000.0);
  EXPECT_EQ(state["STUB2.MODE2.NUMBER"], 500.0);
}

/// @brief Test solving with stub aerosol models
template<class BuilderPolicy>
void test_solve_with_two_stub_aerosol_models(BuilderPolicy builder, double base_relative_tolerance = 5e-5)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1)
                       .AddExternalModelProcesses(aerosol_2)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  // Get a state and set some initial values
  auto state = solver.GetState();

  double fo2_initial = 1.0;
  double baz_mode1_initial = 0.5;
  double fo2_mode2_initial = 0.8;
  double stub2_mode2_fo2_initial = 0.4;
  double stub2_mode2_baz_initial = 0.1;
  double stub2_mode3_baz_initial = 0.2;
  double stub2_mode3_qux_initial = 0.3;
  double temperature = 275.0;
  double fo2_to_baz_rate_constant = 0.01;
  state["FO2"] = fo2_initial;
  state["BAR"] = 2.0;
  state["STUB1.MODE1.QUUX.BAZ"] = baz_mode1_initial;
  state["STUB1.MODE2.CORGE.FO2"] = fo2_mode2_initial;
  state["STUB2.MODE2.CORGE.FO2"] = stub2_mode2_fo2_initial;
  state["STUB2.MODE2.CORGE.BAZ"] = stub2_mode2_baz_initial;
  state["STUB2.MODE3.QUUX.BAZ"] = stub2_mode3_baz_initial;
  state["STUB2.MODE3.QUUX.QUX"] = stub2_mode3_qux_initial;
  state["STUB2.MODE3.CORGE.BAZ"] = 0.2;
  state["STUB2.MODE1.NUMBER"] = 1000.0;
  state["STUB2.MODE2.NUMBER"] = 500.0;

  auto param_it = state.custom_rate_parameter_map_.find("STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
  ASSERT_NE(param_it, state.custom_rate_parameter_map_.end());
  state.custom_rate_parameters_[0][param_it->second] = fo2_to_baz_rate_constant;
  state.conditions_[0].temperature_ = temperature;

  // Calculate the analytical solution to verify the results
  double time_step = 10.0; // seconds
  double stub1_rxn1_delta = fo2_initial * (1.0 - std::exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
  double stub1_rxn2_delta = baz_mode1_initial * (1.0 - std::exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
  double stub2_rxn1_delta = stub2_mode2_fo2_initial * (1.0 - std::exp(-fo2_to_baz_rate_constant * time_step));
  double stub2_rxn2_delta = stub2_mode3_baz_initial * (1.0 - std::exp(-temperature * 0.005 * time_step));
  
  // Solve the system for a single time step
  solver.CalculateRateConstants(state);
  auto results = solver.Solve(time_step, state);

  // Make sure the solver reports success
  EXPECT_EQ(results.state_, micm::SolverState::Converged);

  // Verify that the state variables have been updated
  EXPECT_NEAR(state["FO2"], fo2_initial - stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  EXPECT_EQ(state["BAR"], 2.0);
  EXPECT_NEAR(state["STUB1.MODE1.QUUX.BAZ"], baz_mode1_initial - stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  EXPECT_NEAR(state["STUB1.MODE2.QUUX.BAZ"], stub1_rxn2_delta, base_relative_tolerance * baz_mode1_initial);
  EXPECT_NEAR(state["STUB1.MODE2.CORGE.FO2"], fo2_mode2_initial + stub1_rxn1_delta, base_relative_tolerance * fo2_initial);
  EXPECT_NEAR(state["STUB2.MODE2.CORGE.FO2"], stub2_mode2_fo2_initial - stub2_rxn1_delta, base_relative_tolerance * stub2_mode2_fo2_initial);
  EXPECT_NEAR(state["STUB2.MODE2.CORGE.BAZ"], stub2_mode2_baz_initial + stub2_rxn1_delta, base_relative_tolerance * stub2_mode2_fo2_initial);
  EXPECT_NEAR(state["STUB2.MODE3.QUUX.BAZ"], stub2_mode3_baz_initial - stub2_rxn2_delta, base_relative_tolerance * stub2_mode3_baz_initial);
  EXPECT_NEAR(state["STUB2.MODE3.QUUX.QUX"], stub2_mode3_qux_initial + stub2_rxn2_delta, base_relative_tolerance * stub2_mode3_baz_initial);
  EXPECT_EQ(state["STUB2.MODE1.NUMBER"], 1000.0);
  EXPECT_EQ(state["STUB2.MODE2.NUMBER"], 500.0);
}

/// @brief Test solving with stub aerosol model using 3 grid cells
template<class BuilderPolicy>
void test_solve_with_stub_aerosol_model_1_multi_cell(BuilderPolicy builder, double base_relative_tolerance = 5e-5)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1) // excluding aerosol 2 process for this test
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  const std::size_t num_cells = 3;

  // Get a state for multiple grid cells
  auto state = solver.GetState(num_cells);

  // Set different initial values for each cell
  std::vector<double> fo2_initial = { 1.0, 1.5, 2.0 };
  std::vector<double> baz_mode1_initial = { 0.5, 0.7, 0.9 };
  std::vector<double> fo2_mode2_initial = { 0.8, 1.0, 1.2 };
  
  state["FO2"] = fo2_initial;
  state["BAR"] = std::vector{ 2.0, 2.5, 3.0 };
  state["STUB1.MODE1.QUUX.BAZ"] = baz_mode1_initial;
  state["STUB1.MODE2.CORGE.FO2"] = fo2_mode2_initial;
  state["STUB2.MODE3.CORGE.QUX"] = std::vector{ 0.3, 0.4, 0.5 };
  state["STUB2.MODE3.CORGE.BAZ"] = std::vector{ 0.2, 0.3, 0.4 };
  state["STUB2.MODE1.NUMBER"] = std::vector{ 1000.0, 1100.0, 1200.0 };
  state["STUB2.MODE2.NUMBER"] = std::vector{ 500.0, 550.0, 600.0 };
  
  // Calculate the analytical solution to verify the results for each cell
  double time_step = 10.0; // seconds
  std::vector<double> stub1_rxn1_delta(num_cells);
  std::vector<double> stub1_rxn2_delta(num_cells);
  
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    stub1_rxn1_delta[i] = fo2_initial[i] * (1.0 - std::exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
    stub1_rxn2_delta[i] = baz_mode1_initial[i] * (1.0 - std::exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
  }
  
  // Solve the system for a single time step
  solver.CalculateRateConstants(state);
  auto results = solver.Solve(time_step, state);

  // Make sure the solver reports success
  EXPECT_EQ(results.state_, micm::SolverState::Converged);

  // Verify that the state variables have been updated in all cells
  auto fo2_result = state["FO2"];
  auto bar_result = state["BAR"];
  auto baz_mode1_result = state["STUB1.MODE1.QUUX.BAZ"];
  auto baz_mode2_result = state["STUB1.MODE2.QUUX.BAZ"];
  auto fo2_mode2_result = state["STUB1.MODE2.CORGE.FO2"];
  auto stub2_qux_result = state["STUB2.MODE3.CORGE.QUX"];
  auto stub2_baz_result = state["STUB2.MODE3.CORGE.BAZ"];
  auto stub2_num1_result = state["STUB2.MODE1.NUMBER"];
  auto stub2_num2_result = state["STUB2.MODE2.NUMBER"];
  
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    EXPECT_NEAR(fo2_result[i], fo2_initial[i] - stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    EXPECT_EQ(bar_result[i], (std::vector{ 2.0, 2.5, 3.0 })[i]);
    EXPECT_NEAR(baz_mode1_result[i], baz_mode1_initial[i] - stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    EXPECT_NEAR(baz_mode2_result[i], stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    EXPECT_NEAR(fo2_mode2_result[i], fo2_mode2_initial[i] + stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    EXPECT_EQ(stub2_qux_result[i], (std::vector{ 0.3, 0.4, 0.5 })[i]);
    EXPECT_EQ(stub2_baz_result[i], (std::vector{ 0.2, 0.3, 0.4 })[i]);
    EXPECT_EQ(stub2_num1_result[i], (std::vector{ 1000.0, 1100.0, 1200.0 })[i]);
    EXPECT_EQ(stub2_num2_result[i], (std::vector{ 500.0, 550.0, 600.0 })[i]);
  }
}

/// @brief Test solving with two stub aerosol models using 3 grid cells
template<class BuilderPolicy>
void test_solve_with_two_stub_aerosol_models_multi_cell(BuilderPolicy builder, double base_relative_tolerance = 5e-5)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system with processes that use the aerosol models
  auto solver = builder.SetSystem(system)
                       .AddExternalModelProcesses(aerosol_1)
                       .AddExternalModelProcesses(aerosol_2)
                       .SetIgnoreUnusedSpecies(true)
                       .Build();

  const std::size_t num_cells = 3;

  // Get a state for multiple grid cells
  auto state = solver.GetState(num_cells);

  // Set different initial values for each cell
  std::vector<double> fo2_initial = { 1.0, 1.5, 2.0 };
  std::vector<double> baz_mode1_initial = { 0.5, 0.7, 0.9 };
  std::vector<double> fo2_mode2_initial = { 0.8, 1.0, 1.2 };
  std::vector<double> stub2_mode2_fo2_initial = { 0.4, 0.6, 0.8 };
  std::vector<double> stub2_mode2_baz_initial = { 0.1, 0.15, 0.2 };
  std::vector<double> stub2_mode3_baz_initial = { 0.2, 0.3, 0.4 };
  std::vector<double> stub2_mode3_qux_initial = { 0.3, 0.4, 0.5 };
  std::vector<double> temperature = { 275.0, 285.0, 295.0 };
  std::vector<double> fo2_to_baz_rate_constant = { 0.01, 0.015, 0.02 };
  
  state["FO2"] = fo2_initial;
  state["BAR"] = std::vector{ 2.0, 2.5, 3.0 };
  state["STUB1.MODE1.QUUX.BAZ"] = baz_mode1_initial;
  state["STUB1.MODE2.CORGE.FO2"] = fo2_mode2_initial;
  state["STUB2.MODE2.CORGE.FO2"] = stub2_mode2_fo2_initial;
  state["STUB2.MODE2.CORGE.BAZ"] = stub2_mode2_baz_initial;
  state["STUB2.MODE3.QUUX.BAZ"] = stub2_mode3_baz_initial;
  state["STUB2.MODE3.QUUX.QUX"] = stub2_mode3_qux_initial;
  state["STUB2.MODE3.CORGE.BAZ"] = std::vector{ 0.2, 0.25, 0.3 };
  state["STUB2.MODE1.NUMBER"] = std::vector{ 1000.0, 1100.0, 1200.0 };
  state["STUB2.MODE2.NUMBER"] = std::vector{ 500.0, 550.0, 600.0 };

  auto param_it = state.custom_rate_parameter_map_.find("STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
  ASSERT_NE(param_it, state.custom_rate_parameter_map_.end());
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    state.custom_rate_parameters_[i][param_it->second] = fo2_to_baz_rate_constant[i];
    state.conditions_[i].temperature_ = temperature[i];
  }

  // Calculate the analytical solution to verify the results for each cell
  double time_step = 10.0; // seconds
  std::vector<double> stub1_rxn1_delta(num_cells);
  std::vector<double> stub1_rxn2_delta(num_cells);
  std::vector<double> stub2_rxn1_delta(num_cells);
  std::vector<double> stub2_rxn2_delta(num_cells);
  
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    stub1_rxn1_delta[i] = fo2_initial[i] * (1.0 - std::exp(-STUB1_RATE_CONSTANT_FO2_CORGE * time_step));
    stub1_rxn2_delta[i] = baz_mode1_initial[i] * (1.0 - std::exp(-STUB1_RATE_CONSTANT_BAZ_QUUX * time_step));
    stub2_rxn1_delta[i] = stub2_mode2_fo2_initial[i] * (1.0 - std::exp(-fo2_to_baz_rate_constant[i] * time_step));
    stub2_rxn2_delta[i] = stub2_mode3_baz_initial[i] * (1.0 - std::exp(-temperature[i] * 0.005 * time_step));
  }
  
  // Solve the system for a single time step
  solver.CalculateRateConstants(state);
  auto results = solver.Solve(time_step, state);

  // Make sure the solver reports success
  EXPECT_EQ(results.state_, micm::SolverState::Converged);

  // Verify that the state variables have been updated in all cells
  auto fo2_result = state["FO2"];
  auto bar_result = state["BAR"];
  auto baz_mode1_result = state["STUB1.MODE1.QUUX.BAZ"];
  auto baz_mode2_result = state["STUB1.MODE2.QUUX.BAZ"];
  auto fo2_mode2_result = state["STUB1.MODE2.CORGE.FO2"];
  auto stub2_mode2_fo2_result = state["STUB2.MODE2.CORGE.FO2"];
  auto stub2_mode2_baz_result = state["STUB2.MODE2.CORGE.BAZ"];
  auto stub2_mode3_baz_result = state["STUB2.MODE3.QUUX.BAZ"];
  auto stub2_mode3_qux_result = state["STUB2.MODE3.QUUX.QUX"];
  auto stub2_num1_result = state["STUB2.MODE1.NUMBER"];
  auto stub2_num2_result = state["STUB2.MODE2.NUMBER"];
  
  for (std::size_t i = 0; i < num_cells; ++i)
  {
    EXPECT_NEAR(fo2_result[i], fo2_initial[i] - stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    EXPECT_EQ(bar_result[i], (std::vector{ 2.0, 2.5, 3.0 })[i]);
    EXPECT_NEAR(baz_mode1_result[i], baz_mode1_initial[i] - stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    EXPECT_NEAR(baz_mode2_result[i], stub1_rxn2_delta[i], base_relative_tolerance * baz_mode1_initial[i]);
    EXPECT_NEAR(fo2_mode2_result[i], fo2_mode2_initial[i] + stub1_rxn1_delta[i], base_relative_tolerance * fo2_initial[i]);
    EXPECT_NEAR(stub2_mode2_fo2_result[i], stub2_mode2_fo2_initial[i] - stub2_rxn1_delta[i], base_relative_tolerance * stub2_mode2_fo2_initial[i]);
    EXPECT_NEAR(stub2_mode2_baz_result[i], stub2_mode2_baz_initial[i] + stub2_rxn1_delta[i], base_relative_tolerance * stub2_mode2_fo2_initial[i]);
    EXPECT_NEAR(stub2_mode3_baz_result[i], stub2_mode3_baz_initial[i] - stub2_rxn2_delta[i], base_relative_tolerance * stub2_mode3_baz_initial[i]);
    EXPECT_NEAR(stub2_mode3_qux_result[i], stub2_mode3_qux_initial[i] + stub2_rxn2_delta[i], base_relative_tolerance * stub2_mode3_baz_initial[i]);
    EXPECT_EQ(stub2_num1_result[i], (std::vector{ 1000.0, 1100.0, 1200.0 })[i]);
    EXPECT_EQ(stub2_num2_result[i], (std::vector{ 500.0, 550.0, 600.0 })[i]);
  }
}

