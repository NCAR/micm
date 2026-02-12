// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/test_aerosol_model.cpp
/// @brief Integration test for including an external aerosol model in MICM
///
/// The test uses a stub aerosol model to test the integration interface.
#include <micm/CPU.hpp>

#include <gtest/gtest.h>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <vector>

// some parameters used in the test
constexpr double STUB1_RATE_CONSTANT_FO2_CORGE = 1e-3;
constexpr double STUB1_RATE_CONSTANT_BAZ_QUUX = 2e-3;

// First stubbed aerosol model implementation
//
// This model mimics a single moment two-mode aerosol model
// The first mode contains only the first phase, and the second mode contains both phases
class StubAerosolModel
{
public:
  struct RateConstants
  {
    double fo2_gas_to_mode2_corge; // rate constant for FO2 gas to mode 2 CORGE partitioning
    double baz_mode1_to_mode2_quux; // rate constant for baz mode 1 to baz mode 2 conversion
  };
  StubAerosolModel() = delete;
  StubAerosolModel(const std::string& name, const std::vector<micm::Phase>& phases, const RateConstants& rate_constants) : name_(name), phases_(phases), rate_constants_(rate_constants) {}
  std::tuple<std::size_t, std::size_t> StateSize() const
  {
    EXPECT_EQ(phases_.size(), 2);
    // First mode: first phase only
    // Second mode: both phases
    std::size_t size = 0;
    size += phases_[0].StateSize(); // mode 1
    size += phases_[0].StateSize(); // mode 2, first phase
    size += phases_[1].StateSize(); // mode 2, second phase
    return { size, 0 }; // Return the number of state variables and parameters (0 for this stub model)
  }
  std::vector<std::string> StateVariableNames() const
  {
    std::vector<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    for (const auto& name : phase1_names)
      names.push_back(name_ + ".MODE1." + name);
    for (const auto& name : phase1_names)
      names.push_back(name_ + ".MODE2." + name);
    for (const auto& name : phase2_names)
      names.push_back(name_ + ".MODE2." + name);
    return names;
  }
  std::vector<std::string> StateParameterNames() const
  {
    return {};
  }
  std::string Species(const int mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }

  // We'll assume this model includes gas-aerosol conversion of FO2 to mode 2, and
  // an aerosol-aerosol conversion of baz from mode 1 to mode 2

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    // FO2 gas to mode 2 condensed in CORGE
    auto fo2_gas_index_it = state_indices.find("FO2");
    auto fo2_mode2_index_it = state_indices.find("STUB1.MODE2.CORGE.FO2");
    if (fo2_gas_index_it != state_indices.end() && fo2_mode2_index_it != state_indices.end())
    {
      elements.insert({ fo2_mode2_index_it->second, fo2_gas_index_it->second });
    }
    // baz mode 1 to baz mode 2 in QUUX
    auto baz_mode1_index_it = state_indices.find("STUB1.MODE1.QUUX.BAZ");
    auto baz_mode2_index_it = state_indices.find("STUB1.MODE2.QUUX.BAZ");
    if (baz_mode1_index_it != state_indices.end() && baz_mode2_index_it != state_indices.end())
    {
      elements.insert({ baz_mode2_index_it->second, baz_mode1_index_it->second });
    }
    return elements;
  }
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const DenseMatrixPolicy& state_parameters,
    const DenseMatrixPolicy& state_variables) const
  {
    // We'll store the information needed to calculate the forcing terms in a vector of tuples
    // Each tuple will include: reactant state variable index, product state variable index, and the rate constant
    std::vector<std::tuple<std::size_t, std::size_t, double>> forcing_info;
    auto fo2_gas_index_it = state_variable_indices.find("FO2");
    auto fo2_mode2_index_it = state_variable_indices.find("STUB1.MODE2.CORGE.FO2");
    if (fo2_gas_index_it != state_variable_indices.end() && fo2_mode2_index_it != state_variable_indices.end())
    {
      forcing_info.push_back({ fo2_gas_index_it->second, fo2_mode2_index_it->second, rate_constants_.fo2_gas_to_mode2_corge });
    }
    auto baz_mode1_index_it = state_variable_indices.find("STUB1.MODE1.QUUX.BAZ");
    auto baz_mode2_index_it = state_variable_indices.find("STUB1.MODE2.QUUX.BAZ");
    if (baz_mode1_index_it != state_variable_indices.end() && baz_mode2_index_it != state_variable_indices.end())
    {
      forcing_info.push_back({ baz_mode1_index_it->second, baz_mode2_index_it->second, rate_constants_.baz_mode1_to_mode2_quux });
    }

    // copy-capture the forcing_info vector in the lambda function that will calculate the forcing terms
    return [forcing_info](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing_terms) {
      for (const auto& [reactant_index, product_index, rate_constant] : forcing_info)
      {
        // We'll naively assume the underlying forcing vector is column-major
        // the square-bracket syntax is always [grid_cell][variable_index] regardless of the actual memory layout of the DenseMatrixPolicy
        for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
        {
          // Subtract from reactant
          forcing_terms[i_cell][reactant_index] -= rate_constant * state_variables[i_cell][reactant_index];
          // Add to product
          forcing_terms[i_cell][product_index] += rate_constant * state_variables[i_cell][reactant_index];
        }
      }
    };
  }
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const DenseMatrixPolicy& state_parameters,
    const DenseMatrixPolicy& state_variables,
    const SparseMatrixPolicy& jacobian) const
  {
    // For the simple reactions in this stub model, there will be two Jacobian element updates for each reaction: one for the reactant and one for the product, with values equal to the rate constant (negative for the reactant, positive for the product)
    // We'll store the information needed to calculate the Jacobian terms in a vector of tuples, similar to the forcing function
    // with jacobian index, and the partial derivative value (+/-rate constant)
    std::vector<std::tuple<std::size_t, double>> jacobian_info;
    for(std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
    {
      auto fo2_gas_index_it = state_variable_indices.find("FO2");
      auto fo2_mode2_index_it = state_variable_indices.find("STUB1.MODE2.CORGE.FO2");
      if (fo2_gas_index_it != state_variable_indices.end() && fo2_mode2_index_it != state_variable_indices.end())
      {
        jacobian_info.push_back({ jacobian.VectorIndex(i_cell, fo2_gas_index_it->second, fo2_gas_index_it->second), -rate_constants_.fo2_gas_to_mode2_corge }); // reactant partial derivative
        jacobian_info.push_back({ jacobian.VectorIndex(i_cell, fo2_mode2_index_it->second, fo2_gas_index_it->second), rate_constants_.fo2_gas_to_mode2_corge }); // product partial derivative
      }
      auto baz_mode1_index_it = state_variable_indices.find("STUB1.MODE1.QUUX.BAZ");
      auto baz_mode2_index_it = state_variable_indices.find("STUB1.MODE2.QUUX.BAZ");
      if (baz_mode1_index_it != state_variable_indices.end() && baz_mode2_index_it != state_variable_indices.end())
      {
        jacobian_info.push_back({ jacobian.VectorIndex(i_cell, baz_mode1_index_it->second, baz_mode1_index_it->second), -rate_constants_.baz_mode1_to_mode2_quux }); // reactant partial derivative
        jacobian_info.push_back({ jacobian.VectorIndex(i_cell, baz_mode2_index_it->second, baz_mode1_index_it->second), rate_constants_.baz_mode1_to_mode2_quux }); // product partial derivative
      } 
    }

    // copy-capture the jacobian_info vector in the lambda function that will calculate the Jacobian terms
    return [jacobian_info](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) {
      for (const auto& [jacobian_index, value] : jacobian_info)
      {
        jacobian.AsVector()[jacobian_index] += value;
      }
    };
  }

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
  RateConstants rate_constants_;
};

// Second stubbed aerosol model implementation
//
// This model mimics a two-moment three-mode aerosol model
// The first mode contains only the first phase, the second mode contains the second phase,
// and the third mode contains both phases. For two-moment schemes, particle number concentration
// is also included for each mode.
class AnotherStubAerosolModel
{
public:
  AnotherStubAerosolModel() = delete;
  AnotherStubAerosolModel(const std::string& name, const std::vector<micm::Phase>& phases) : name_(name), phases_(phases) {}
  std::tuple<std::size_t, std::size_t> StateSize() const
  {
    EXPECT_EQ(phases_.size(), 2);
    // First mode: first phase only
    // Second mode: second phase only
    // Third mode: both phases
    std::size_t size = 0;
    size += 1; // mode 1 number concentration
    size += phases_[0].StateSize(); // mode 1 species
    size += 1; // mode 2 number concentration
    size += phases_[1].StateSize(); // mode 2 species
    size += 1; // mode 3 number concentration
    size += phases_[0].StateSize(); // mode 3, first phase species
    size += phases_[1].StateSize(); // mode 3, second phase species
    return { size, 0 }; // Return the number of state variables and parameters (0 for this stub model for now)
  }
  std::vector<std::string> StateVariableNames() const
  {
    std::vector<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    names.push_back(name_ + ".MODE1.NUMBER"); // number concentration for mode 1
    for (const auto& name : phase1_names)
    names.push_back(name_ + ".MODE1." + name);
    names.push_back(name_ + ".MODE2.NUMBER"); // number concentration for mode 2
    for (const auto& name : phase2_names)
    names.push_back(name_ + ".MODE2." + name);
    names.push_back(name_ + ".MODE3.NUMBER"); // number concentration for mode 3
    for (const auto& name : phase1_names)
    names.push_back(name_ + ".MODE3." + name);
    for (const auto& name : phase2_names)
    names.push_back(name_ + ".MODE3." + name);
    return names;
  }
  std::vector<std::string> StateParameterNames() const
  {
    return {};
  }
  std::string Species(const int mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }
  std::string Number(const int mode) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + ".NUMBER";
  }
  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    // For this stub model, we'll assume there are no processes and therefore the Jacobian is always zero, so we can return an empty set
    return {};
  }
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const DenseMatrixPolicy& state_parameters,
    const DenseMatrixPolicy& state_variables) const
  {
    // For this stub model, we'll assume there are no processes and therefore the forcing function always returns zero, so we can return a lambda that does nothing
    return [](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing_terms) {
      // Do nothing
    };
  }
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const DenseMatrixPolicy& state_parameters,
    const DenseMatrixPolicy& state_variables,
    const SparseMatrixPolicy& jacobian) const
  {
    // For this stub model, we'll assume there are no processes and therefore the Jacobian is always zero, so we can return a lambda that does nothing
    return [](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) {
      // Do nothing
    };
  }

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
};

std::tuple<micm::System, StubAerosolModel, AnotherStubAerosolModel, std::map<std::string, micm::Phase>> CreateSystemWithStubAerosolModels()
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

TEST(AerosolModelIntegration, CanIntegrateWithStubAerosolModel)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  micm::Solver solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(system)
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

TEST(AerosolModelIntegration, CanUpdateStateWithStubAerosolModel)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  micm::Solver solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(system)
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

TEST(AerosolModelIntegration, CanUpdateMultiCellStateWithStubAerosolModel)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  micm::Solver solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(system)
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

TEST(AerosolModelIntegration, CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  micm::Solver solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(system)
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
  auto forcing_function_1 = aerosol_1.ForcingFunction(state.custom_rate_parameter_map_, state.variable_map_, state.custom_rate_parameters_, state.variables_);
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

TEST(AerosolModelIntegration, CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
  auto [system, aerosol_1, aerosol_2, phases] = CreateSystemWithStubAerosolModels();

  // Create a solver for the system (without processes for simplicity)
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  micm::Solver solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(system)
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
  auto jacobian_function_1 = aerosol_1.JacobianFunction(state.custom_rate_parameter_map_, state.variable_map_, state.custom_rate_parameters_, state.variables_, jacobian_1);
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
  EXPECT_DOUBLE_EQ(jacobian_1[0][fo2_gas_index][fo2_gas_index], expected_fo2_gas_partial);
  EXPECT_DOUBLE_EQ(jacobian_1[0][fo2_mode2_index][fo2_gas_index], expected_fo2_mode2_partial);

  // For the baz mode 1 to mode 2 QUUX conversion, we expect two non-zero Jacobian elements: a negative value equal to the rate constant in the column corresponding to baz mode 1 QUUX and row corresponding to baz mode 1 QUUX (representing the partial derivative of the baz loss with respect to baz), and a positive value equal to the rate constant in the column corresponding to baz mode 1 QUUX and row corresponding to baz mode 2 QUUX (representing the partial derivative of the baz gain with respect to baz)
  auto baz_mode1_index_it = state.variable_map_.find("STUB1.MODE1.QUUX.BAZ");
  auto baz_mode2_index_it = state.variable_map_.find("STUB1.MODE2.QUUX.BAZ");
  ASSERT_NE(baz_mode1_index_it, state.variable_map_.end());
  ASSERT_NE(baz_mode2_index_it, state.variable_map_.end());
  std::size_t baz_mode1_index = baz_mode1_index_it->second;
  std::size_t baz_mode2_index = baz_mode2_index_it->second;
  double expected_baz_mode1_partial = -STUB1_RATE_CONSTANT_BAZ_QUUX;
  double expected_baz_mode2_partial = STUB1_RATE_CONSTANT_BAZ_QUUX;
  EXPECT_DOUBLE_EQ(jacobian_1[0][baz_mode1_index][baz_mode1_index], expected_baz_mode1_partial);
  EXPECT_DOUBLE_EQ(jacobian_1[0][baz_mode2_index][baz_mode1_index], expected_baz_mode2_partial);
}