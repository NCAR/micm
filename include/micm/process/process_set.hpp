// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <vector>

#ifdef USE_CUDA
#include <micm/process/process_set_cuda.hpp>
#endif

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
    std::vector<std::size_t> jacobian_flat_ids_;

   public:
    /// @brief Default constructor
    ProcessSet() = default;

    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param state Solver state
    ProcessSet(const std::vector<Process>& processes, const State& state);

    /// @brief Return the full set of non-zero Jacobian elements for the set of processes
    /// @return Jacobian elements as a set of index pairs
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;

    /// @brief Sets the indicies for each non-zero Jacobian element in the underlying vector
    /// @param matrix The sparse matrix used for the Jacobian
    void SetJacobianFlatIds(const SparseMatrix<double>& matrix);

    /// @brief Add forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    void AddForcingTerms(
        const Matrix<double>& rate_constants,
        const Matrix<double>& state_variables,
        Matrix<double>& forcing) const;

    /// @brief Add Jacobian terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param jacobian Jacobian matrix for the system (grid cell, dependent variable, independent variable)
    void AddJacobianTerms(
        const Matrix<double>& rate_constants,
        const Matrix<double>& state_variables,
        SparseMatrix<double>& jacobian) const;
    
    #ifdef USE_CUDA
    void CudaAddForcingTerms(
        const Matrix<double>& rate_constants, 
        const Matrix<double>& state_variables, 
        Matrix<double>& forcing);
    #endif
  
  };

  inline ProcessSet::ProcessSet(const std::vector<Process>& processes, const State& state)
      : number_of_reactants_(),
        reactant_ids_(),
        number_of_products_(),
        product_ids_(),
        yields_()
  {
    for (auto& process : processes)
    {
      number_of_reactants_.push_back(process.reactants_.size());
      number_of_products_.push_back(process.products_.size());
      for (auto& reactant : process.reactants_)
      {
        reactant_ids_.push_back(state.variable_map_.at(reactant.name_));
      }
      for (auto& product : process.products_)
      {
        product_ids_.push_back(state.variable_map_.at(product.first.name_));
        yields_.push_back(product.second);
      }
    }
  };

  inline std::set<std::pair<std::size_t, std::size_t>> ProcessSet::NonZeroJacobianElements() const
  {
    std::set<std::pair<std::size_t, std::size_t>> ids;
    auto react_id = reactant_ids_.begin();
    auto prod_id = product_ids_.begin();
    for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      for (std::size_t i_ind = 0; i_ind < number_of_reactants_[i_rxn]; ++i_ind)
      {
        for (std::size_t i_dep = 0; i_dep < number_of_reactants_[i_rxn]; ++i_dep)
        {
          ids.insert(std::make_pair(react_id[i_dep], react_id[i_ind]));
        }
        for (std::size_t i_dep = 0; i_dep < number_of_products_[i_rxn]; ++i_dep)
        {
          ids.insert(std::make_pair(prod_id[i_dep], react_id[i_ind]));
        }
      }
      react_id += number_of_reactants_[i_rxn];
      prod_id += number_of_products_[i_rxn];
    }
    return ids;
  }

  inline void ProcessSet::SetJacobianFlatIds(const SparseMatrix<double>& matrix)
  {
    jacobian_flat_ids_.clear();
    auto react_id = reactant_ids_.begin();
    auto prod_id = product_ids_.begin();
    for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      for (std::size_t i_ind = 0; i_ind < number_of_reactants_[i_rxn]; ++i_ind)
      {
        for (std::size_t i_dep = 0; i_dep < number_of_reactants_[i_rxn]; ++i_dep)
        {
          jacobian_flat_ids_.push_back(matrix.VectorIndex(0, react_id[i_dep], react_id[i_ind]));
        }
        for (std::size_t i_dep = 0; i_dep < number_of_products_[i_rxn]; ++i_dep)
        {
          jacobian_flat_ids_.push_back(matrix.VectorIndex(0, prod_id[i_dep], react_id[i_ind]));
        }
      }
      react_id += number_of_reactants_[i_rxn];
      prod_id += number_of_products_[i_rxn];
    }
  }

  inline void ProcessSet::AddForcingTerms(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      Matrix<double>& forcing) const
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

  inline void ProcessSet::AddJacobianTerms(
      const Matrix<double>& rate_constants,
      const Matrix<double>& state_variables,
      SparseMatrix<double>& jacobian) const
  {
    auto cell_jacobian = jacobian.AsVector().begin();
    // loop over grid cells
    for (std::size_t i_cell = 0; i_cell < state_variables.size(); ++i_cell)
    {
      auto cell_rate_constants = rate_constants[i_cell];
      auto cell_state = state_variables[i_cell];
      auto react_id = reactant_ids_.begin();
      auto yield = yields_.begin();
      auto flat_id = jacobian_flat_ids_.begin();
      for (std::size_t i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
      {
        for (std::size_t i_ind = 0; i_ind < number_of_reactants_[i_rxn]; ++i_ind)
        {
          double d_rate_d_ind = cell_rate_constants[i_rxn];
          for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
          {
            if (i_react == i_ind)
              continue;
            d_rate_d_ind *= cell_state[react_id[i_react]];
          }
          for (std::size_t i_dep = 0; i_dep < number_of_reactants_[i_rxn]; ++i_dep)
            cell_jacobian[*(flat_id++)] -= d_rate_d_ind;
          for (std::size_t i_dep = 0; i_dep < number_of_products_[i_rxn]; ++i_dep)
            cell_jacobian[*(flat_id++)] += yield[i_dep] * d_rate_d_ind;
        }
        react_id += number_of_reactants_[i_rxn];
        yield += number_of_products_[i_rxn];
      }
      cell_jacobian += jacobian.FlatBlockSize();
    }
  }

    #ifdef USE_CUDA
    inline void ProcessSet::CudaAddForcingTerms(
        const Matrix<double>& rate_constants, 
        const Matrix<double>& state_variables, 
        Matrix<double>& forcing) {
          AddForcingTerms_kernelSetup(
            rate_constants, 
            state_variables, 
            forcing,
            number_of_reactants_,
            reactant_ids_,
            number_of_products_,
            product_ids_,
            yields_,
            jacobian_flat_ids_
          );
        }
    #endif
}  // namespace micm