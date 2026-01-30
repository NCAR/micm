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
    EXPECT_TRUE(phases_.size() == 2);
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
    EXPECT_TRUE(phases_.size() == 2);
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
  AnotherStubAerosolModel(std::string name, std::vector<micm::Phase> phases) : name_(name), phases_(phases) {}
  std::size_t StateSize() const
  {
    EXPECT_TRUE(phases_.size() == 2);
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
    EXPECT_TRUE(phases_.size() == 2);
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

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
};

TEST(AerosolModelIntegration, CanIntegrateWithStubAerosolModel)
{
  // Create a simple chemical system
  auto foo = micm::Species("FO2"); // species that can partition to condensed phase
  auto bar = micm::Species("BAR"); // gas-phase only species
  auto baz = micm::Species("BAZ"); // condensed-phase only species
  auto qux = micm::Species("QUX"); // condensed-phase only species

  auto gas   = micm::Phase("GAS", std::vector<micm::PhaseSpecies>({ foo, bar }));        // gas phase
  auto quux  = micm::Phase("QUUX", std::vector<micm::PhaseSpecies>({ baz, qux }));       // condensed aerosol or cloud phase
  auto corge = micm::Phase("CORGE", std::vector<micm::PhaseSpecies>({ foo, baz, qux })); // another condensed aerosol or cloud phase

  // Create instances of each stub aerosol model
  auto aerosol_1 = StubAerosolModel("STUB1", std::vector<micm::Phase>({ quux, corge }));
  auto aerosol_2 = AnotherStubAerosolModel("STUB2", std::vector<micm::Phase>({ quux, corge }));

  // Create a system containing the gas phase and both aerosol models
  auto system = micm::System({
    .gas_phase_ = gas,
    .external_models_ = { aerosol_1, aerosol_2 }
  });

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
      gas.UniqueNames().size() + aerosol_1.UniqueNames().size() + aerosol_2.UniqueNames().size());

  // Assemble the full list of expected variable names
  std::vector<std::string> expected_names;
  auto gas_names = gas.SpeciesNames();
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