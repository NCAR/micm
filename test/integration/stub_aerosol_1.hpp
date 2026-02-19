// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_1.hpp
/// @brief Stub aerosol model for integration testing
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
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
  std::set<std::string> StateVariableNames() const
  {
    std::set<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    for (const auto& name : phase1_names)
      names.insert(name_ + ".MODE1." + name);
    for (const auto& name : phase1_names)
      names.insert(name_ + ".MODE2." + name);
    for (const auto& name : phase2_names)
      names.insert(name_ + ".MODE2." + name);
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

  // We'll pretend all the species in the aerosol model are involved in processes
  std::set<std::string> SpeciesUsed() const
  {
    return StateVariableNames();
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

  // We have no parameters for this stub model
  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
  {
    // No parameters to update in this stub model
    return [](const std::vector<micm::Conditions>& conditions, DenseMatrixPolicy& state_parameters) {
      // Do nothing
    };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
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
    const SparseMatrixPolicy& jacobian) const
  {
    // For this simple implementation, we'll use the dependent and independent variable indices with the square-bracket syntax
    // of the jacobian matrix. In a real implementation, we should want to get the underlying vector indices of the jacobian
    // elements, and iterate over blocks in the block diagonal sparse matrix in the most efficient way for the specific SparseMatrixPolicy.
    std::vector<std::tuple<std::size_t, std::size_t, double>> jacobian_info; // (dependent id, independent id, value)
    auto fo2_gas_index_it = state_variable_indices.find("FO2");
    auto fo2_mode2_index_it = state_variable_indices.find("STUB1.MODE2.CORGE.FO2");
    if (fo2_gas_index_it != state_variable_indices.end() && fo2_mode2_index_it != state_variable_indices.end())
    {
      jacobian_info.push_back({ fo2_gas_index_it->second, fo2_gas_index_it->second, -rate_constants_.fo2_gas_to_mode2_corge }); // reactant partial derivative
      jacobian_info.push_back({ fo2_mode2_index_it->second, fo2_gas_index_it->second, rate_constants_.fo2_gas_to_mode2_corge }); // product partial derivative
    }
    auto baz_mode1_index_it = state_variable_indices.find("STUB1.MODE1.QUUX.BAZ");
    auto baz_mode2_index_it = state_variable_indices.find("STUB1.MODE2.QUUX.BAZ");
    if (baz_mode1_index_it != state_variable_indices.end() && baz_mode2_index_it != state_variable_indices.end())
    {
      jacobian_info.push_back({ baz_mode1_index_it->second, baz_mode1_index_it->second, -rate_constants_.baz_mode1_to_mode2_quux }); // reactant partial derivative
      jacobian_info.push_back({ baz_mode2_index_it->second, baz_mode1_index_it->second, rate_constants_.baz_mode1_to_mode2_quux }); // product partial derivative
    } 
  
    // copy-capture the jacobian_info vector in the lambda function that will calculate the Jacobian terms
    return [jacobian_info](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) {
      for (std::size_t i_block = 0; i_block < jacobian.NumberOfBlocks(); ++i_block)
      {
        for (const auto& [dependent_id, independent_id, value] : jacobian_info)
        {
          jacobian[i_block][dependent_id][independent_id] -= value;
        }
      }
    };
  }

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
  RateConstants rate_constants_;
};
