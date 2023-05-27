// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <vector>

namespace micm
{

  /// @brief Solver function calculators for a collection of processes
  class ProcessSet
  {
    std::vector<std::size_t> number_of_reactants_;
    std::vector<std::size_t> reactant_ids_;
    std::vector<std::size_t> number_of_products_;
    std::vector<std::size_t> product_ids_;
    std::vector<double> yields_;

   public:
    /// @brief Default constructor
    ProcessSet() = default;

    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calcualtor for
    /// @param state Solver state
    ProcessSet(const std::vector<Process>& processes, const State& state);

    /// @brief Add forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    void
    AddForcingTerms(const Matrix<double>& rate_constants, const Matrix<double>& state_variables, Matrix<double>& forcing);
  };

  inline ProcessSet::ProcessSet(const std::vector<Process>& processes, const State& state)
      : number_of_reactants_(
            [&]() -> std::vector<std::size_t>
            {
              std::vector<std::size_t> num{};
              for (auto& process : processes)
              {
                num.push_back(process.reactants_.size());
              }
              return num;
            }()),
        reactant_ids_(
            [&]() -> std::vector<std::size_t>
            {
              std::vector<std::size_t> ids{};
              for (auto& process : processes)
              {
                for (auto& reactant : process.reactants_)
                {
                  ids.push_back(state.variable_map_.at(reactant.name_));
                }
              }
              return ids;
            }()),
        number_of_products_(
            [&]() -> std::vector<std::size_t>
            {
              std::vector<std::size_t> num{};
              for (auto& process : processes)
              {
                num.push_back(process.products_.size());
              }
              return num;
            }()),
        product_ids_(
            [&]() -> std::vector<std::size_t>
            {
              std::vector<std::size_t> ids{};
              for (auto& process : processes)
              {
                for (auto& product : process.products_)
                {
                  ids.push_back(state.variable_map_.at(product.first.name_));
                }
              }
              return ids;
            }()),
        yields_(
            [&]() -> std::vector<double>
            {
              std::vector<double> yields{};
              for (auto& process : processes)
              {
                for (auto& product : process.products_)
                {
                  yields.push_back(product.second);
                }
              }
              return yields;
            }()){};

  inline void ProcessSet::AddForcingTerms(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      Matrix<double>& forcing)
  {
    // loop over grid cells
    for (std::size_t i_cell = 0; i_cell < state_variables.size(); ++i_cell)
    {
      auto cell_rate_constants = rate_constants[i_cell];
      auto cell_state = state_variables[i_cell];
      auto cell_forcing = forcing[i_cell];
      auto react_id = reactant_ids_.begin();
      auto prod_id = product_ids_.begin();
      auto yield = yields_.begin();
      for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
      {
        double rate = cell_rate_constants[i_rxn];
        for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
          rate *= cell_state[react_id[i_react]];
        for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
          cell_forcing[react_id[i_react]] -= rate;
        for (std::size_t i_prod = 0; i_prod < number_of_products_[i_rxn]; ++i_prod)
          cell_forcing[prod_id[i_prod]] += yield[i_prod] * rate;
        react_id += number_of_reactants_[i_rxn];
        prod_id += number_of_products_[i_rxn];
        yield += number_of_products_[i_rxn];
      }
    }
  };

}  // namespace micm