// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_2.hpp
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
    return { size, 2 }; // Return the number of state variables and parameters (2 for this stub model)
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
    std::set<std::string> names;

    // add parameters for the two rate constants used in the model processes
    names.insert(name_ + ".PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
    names.insert(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
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

  // We'll assume all the species in the aerosol model are involved in processes
  std::set<std::string> SpeciesUsed() const
  {
    return StateVariableNames();
  }

  // We'll assume this model includes two processes:
  // 1) conversion of FO2 to BAZ in mode 2 CORGE phase
  // 2) conversion of BAZ to QUX in mode 3 QUUX phase

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    // FO2 to BAZ in mode 2 CORGE
    auto fo2_mode2_index_it = state_indices.find("STUB2.MODE2.CORGE.FO2");
    auto baz_mode2_index_it = state_indices.find("STUB2.MODE2.CORGE.BAZ");
    if (fo2_mode2_index_it != state_indices.end() && baz_mode2_index_it != state_indices.end())
    {
      elements.insert({ baz_mode2_index_it->second, fo2_mode2_index_it->second });
    }
    // BAZ to QUX in mode 3 QUUX
    auto baz_mode3_index_it = state_indices.find("STUB2.MODE3.QUUX.BAZ");
    auto qux_mode3_index_it = state_indices.find("STUB2.MODE3.QUUX.QUX");
    if (baz_mode3_index_it != state_indices.end() && qux_mode3_index_it != state_indices.end())
    {
      elements.insert({ qux_mode3_index_it->second, baz_mode3_index_it->second });
    }
    return elements;
  }

  // We have parameters for this stub model, one of which we will update based on temperature
  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
  {
    // Create a function that updates the BAZ-to-QUX rate constant based on temperature for each grid cell
    auto baz_to_qux_param_it = state_parameter_indices.find(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
    EXPECT_NE(baz_to_qux_param_it, state_parameter_indices.end());
    std::size_t baz_to_qux_param_index = baz_to_qux_param_it->second;
    return [baz_to_qux_param_index](const std::vector<micm::Conditions>& conditions, DenseMatrixPolicy& state_parameters) {
      for(std::size_t cell = 0; cell < conditions.size(); ++cell)
      {
        // Simple linear dependence on temperature for this stub model
        double temperature = conditions[cell].temperature_;
        double rate_constant = 0.005 * temperature; // simple function for testing
        state_parameters[cell][baz_to_qux_param_index] = rate_constant;
      }
    };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
    const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
    const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
  {
    // We'll store the information needed to calculate the forcing terms in a vector of tuples
    // Each tuple will include: reactant state variable index, product state variable index, and the index of the rate constant parameter
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> forcing_info;
    auto fo2_mode2_index_it = state_variable_indices.find("STUB2.MODE2.CORGE.FO2");
    auto baz_mode2_index_it = state_variable_indices.find("STUB2.MODE2.CORGE.BAZ");
    auto fo2_to_baz_param_it = state_parameter_indices.find(name_ + ".PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
    if (fo2_mode2_index_it != state_variable_indices.end() &&
        baz_mode2_index_it != state_variable_indices.end() &&
        fo2_to_baz_param_it != state_parameter_indices.end())
    {
      forcing_info.push_back({ fo2_mode2_index_it->second, baz_mode2_index_it->second, fo2_to_baz_param_it->second });
    }
    auto baz_mode3_index_it = state_variable_indices.find("STUB2.MODE3.QUUX.BAZ");
    auto qux_mode3_index_it = state_variable_indices.find("STUB2.MODE3.QUUX.QUX");
    auto baz_to_qux_param_it = state_parameter_indices.find(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
    if (baz_mode3_index_it != state_variable_indices.end() &&
        qux_mode3_index_it != state_variable_indices.end() &&
        baz_to_qux_param_it != state_parameter_indices.end())
    {
      forcing_info.push_back({ baz_mode3_index_it->second, qux_mode3_index_it->second, baz_to_qux_param_it->second });
    }

    // copy capture the forcing_info vector in the lambda function that will calculate the forcing terms
    return [forcing_info](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing_terms) {
      // We'll naively assume the underlying forcing vector is column-major
      // the square-bracket syntax is always [grid_cell][variable_index] regardless of the actual memory layout of the DenseMatrixPolicy
      for (const auto& [reactant_index, product_index, rate_param_index] : forcing_info)
      {
        for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
        {
          double rate_constant = state_parameters[i_cell][rate_param_index];
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
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, double>> jacobian_info; // (dependent id, independent id, rate param id, multiplier)
    auto fo2_mode2_index_it = state_variable_indices.find("STUB2.MODE2.CORGE.FO2");
    auto baz_mode2_index_it = state_variable_indices.find("STUB2.MODE2.CORGE.BAZ");
    auto fo2_to_baz_param_it = state_parameter_indices.find(name_ + ".PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
    if (fo2_mode2_index_it != state_variable_indices.end() &&
        baz_mode2_index_it != state_variable_indices.end() &&
        fo2_to_baz_param_it != state_parameter_indices.end())
    {
      jacobian_info.push_back({ fo2_mode2_index_it->second, fo2_mode2_index_it->second, fo2_to_baz_param_it->second, -1.0 }); // reactant partial derivative
      jacobian_info.push_back({ baz_mode2_index_it->second, fo2_mode2_index_it->second, fo2_to_baz_param_it->second, 1.0 }); // product partial derivative
    }
    auto baz_mode3_index_it = state_variable_indices.find("STUB2.MODE3.QUUX.BAZ");
    auto qux_mode3_index_it = state_variable_indices.find("STUB2.MODE3.QUUX.QUX");
    auto baz_to_qux_param_it = state_parameter_indices.find(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
    if (baz_mode3_index_it != state_variable_indices.end() &&
        qux_mode3_index_it != state_variable_indices.end() &&
        baz_to_qux_param_it != state_parameter_indices.end())
    {
      jacobian_info.push_back({ baz_mode3_index_it->second, baz_mode3_index_it->second, baz_to_qux_param_it->second, -1.0 }); // reactant partial derivative
      jacobian_info.push_back({ qux_mode3_index_it->second, baz_mode3_index_it->second, baz_to_qux_param_it->second, 1.0 }); // product partial derivative
    }

    // copy-capture the jacobian_info vector in the lambda function that will calculate the Jacobian terms
    return [jacobian_info](const DenseMatrixPolicy& state_parameters, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) {
      for (std::size_t i_block = 0; i_block < jacobian.NumberOfBlocks(); ++i_block)
      {
        for (const auto& [dependent_id, independent_id, rate_param_id, multiplier] : jacobian_info)
        {
          double rate_constant = state_parameters[i_block][rate_param_id];
          jacobian[i_block][dependent_id][independent_id] -= multiplier * rate_constant;
        }
      }
    };
  }

private:
  std::string name_;
  std::vector<micm::Phase> phases_;
};
