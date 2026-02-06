// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/test_aerosol_model.cpp
/// @brief Integration test for including an external aerosol model in MICM
///
/// The test uses a stub aerosol model to test the integration interface.
#include <micm/CPU.hpp>


#include <gtest/gtest.h>

// First stubbed aerosol model implementation
//
// This model mimics a single moment two-mode aerosol model
// The first mode contains only the first phase, and the second mode contains both phases
class StubAerosolModel
{
public:
  StubAerosolModel() = delete;
  StubAerosolModel(const std::string& name, const std::vector<micm::Phase>& phases) : name_(name), phases_(phases) {}
  std::size_t StateSize() const
  {
    EXPECT_EQ(phases_.size(), 2);
    // First mode: first phase only
    // Second mode: both phases
    std::size_t size = 0;
    size += phases_[0].StateSize(); // mode 1
    size += phases_[0].StateSize(); // mode 2, first phase
    size += phases_[1].StateSize(); // mode 2, second phase
    return size;
  }
  std::vector<std::string> UniqueNames() const
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
  std::string Species(const int mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
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
  std::size_t StateSize() const
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
    return size;
  }
  std::vector<std::string> UniqueNames() const
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
  std::string Species(const int mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }
  std::string Number(const int mode) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + ".NUMBER";
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
  auto aerosol_1 = StubAerosolModel("STUB1", std::vector<micm::Phase>({ quux, corge }));
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
      system.gas_phase_.UniqueNames().size() + aerosol_1.UniqueNames().size() + aerosol_2.UniqueNames().size());

  // Assemble the full list of expected variable names
  std::vector<std::string> expected_names;
  auto gas_names = system.gas_phase_.SpeciesNames();
  expected_names.insert(expected_names.end(), gas_names.begin(), gas_names.end());
  auto aerosol1_names = aerosol_1.UniqueNames();
  expected_names.insert(expected_names.end(), aerosol1_names.begin(), aerosol1_names.end());
  auto aerosol2_names = aerosol_2.UniqueNames();
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
  auto& corge = phases["CORGE"];
  auto qux = micm::Species("QUX");
  state[aerosol_2.Species(2, corge, qux)] = 0.42;
  // Set some condensed-phase species in the second aerosol model by index
  auto corge_baz_index = state.variable_map_.find("STUB2.MODE3.CORGE.BAZ")->second;
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
  auto& corge = phases["CORGE"];
  auto qux = micm::Species("QUX");
  state[aerosol_2.Species(2, corge, qux)] = std::vector{ 0.42, 0.43, 0.44 };

  // Set some condensed-phase species in the second aerosol model by index
  auto corge_baz_index = state.variable_map_.find("STUB2.MODE3.CORGE.BAZ")->second;
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