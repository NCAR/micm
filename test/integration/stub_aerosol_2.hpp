// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_1.hpp
/// @brief Stub aerosol model for integration testing
#pragma once

#include <gtest/gtest.h>

#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

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
  std::set<std::string> StateVariableNames() const
  {
    std::set<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    names.insert(name_ + ".MODE1.NUMBER"); // number concentration for mode 1
    for (const auto& name : phase1_names)
    names.insert(name_ + ".MODE1." + name);
    names.insert(name_ + ".MODE2.NUMBER"); // number concentration for mode 2
    for (const auto& name : phase2_names)
    names.insert(name_ + ".MODE2." + name);
    names.insert(name_ + ".MODE3.NUMBER"); // number concentration for mode 3
    for (const auto& name : phase1_names)
    names.insert(name_ + ".MODE3." + name);
    for (const auto& name : phase2_names)
    names.insert(name_ + ".MODE3." + name);
    return names;
  }
  std::set<std::string> StateParameterNames() const
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

  // We'll assume all the species in the aerosol model are involved in processes
  std::set<std::string> SpeciesUsed() const
  {
    return StateVariableNames();
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    // For this stub model, we'll assume there are no processes and therefore the Jacobian is always zero, so we can return an empty set
    return {};
  }
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
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