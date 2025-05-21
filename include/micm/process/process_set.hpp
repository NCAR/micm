// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/error.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cassert>
#include <vector>

enum class MicmProcessSetErrc
{
  ReactantDoesNotExist = MICM_PROCESS_SET_ERROR_CODE_REACTANT_DOES_NOT_EXIST,
  ProductDoesNotExist = MICM_PROCESS_SET_ERROR_CODE_PRODUCT_DOES_NOT_EXIST
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmProcessSetErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmProcessSetErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_PROCESS_SET;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmProcessSetErrc>(ev))
      {
        case MicmProcessSetErrc::ReactantDoesNotExist: return "Reactant does not exist";
        case MicmProcessSetErrc::ProductDoesNotExist: return "Product does not exist";
        default: return "Unknown error";
      }
    }
  };

  const MicmProcessSetErrorCategory MICM_PROCESS_SET_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmProcessSetErrc e)
{
  return { static_cast<int>(e), MICM_PROCESS_SET_ERROR };
}

namespace micm
{

  /// @brief Solver function calculators for a collection of processes
  class ProcessSet
  {
   protected:
    /// @brief Process information for use in setting Jacobian elements
    struct ProcessInfo
    {
      std::size_t process_id_;
      std::size_t independent_id_;
      std::size_t number_of_dependent_reactants_;
      std::size_t number_of_products_;
    };

    std::vector<std::size_t> number_of_reactants_;
    std::vector<std::size_t> reactant_ids_;
    std::vector<std::size_t> number_of_products_;
    std::vector<std::size_t> product_ids_;
    std::vector<double> yields_;
    std::vector<ProcessInfo> jacobian_process_info_;
    std::vector<std::size_t> jacobian_reactant_ids_;
    std::vector<std::size_t> jacobian_product_ids_;
    std::vector<double> jacobian_yields_;
    std::vector<std::size_t> jacobian_flat_ids_;

   public:
    /// @brief Default constructor
    ProcessSet() = default;

    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    ProcessSet(const std::vector<Process>& processes, const std::map<std::string, std::size_t>& variable_map);

    virtual ~ProcessSet() = default;

    /// @brief Return the full set of non-zero Jacobian elements for the set of processes
    /// @return Jacobian elements as a set of index pairs
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;

    /// @brief Sets the indicies for each non-zero Jacobian element in the underlying vector.
    ///        Also sets combination of process ids and reactant ids to allow setting of
    ///        jacobian elements in order by column.
    /// @param matrix The sparse matrix used for the Jacobian
    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);

    /// @brief Add forcing terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    template<typename DenseMatrixPolicy>
      requires(!VectorizableDense<DenseMatrixPolicy>)
    void AddForcingTerms(
        const DenseMatrixPolicy& rate_constants,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const;
    template<typename DenseMatrixPolicy>
      requires VectorizableDense<DenseMatrixPolicy>
    void AddForcingTerms(
        const DenseMatrixPolicy& rate_constants,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const;

    /// @brief Subtract Jacobian terms for the set of processes for the current conditions
    /// @param rate_constants Current values for the process rate constants (grid cell, process)
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param jacobian Jacobian matrix for the system (grid cell, dependent variable, independent variable)
    template<class DenseMatrixPolicy, class SparseMatrixPolicy>
      requires(!VectorizableDense<DenseMatrixPolicy> || !VectorizableSparse<SparseMatrixPolicy>)
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& rate_constants,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const;
    // template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
    template<class DenseMatrixPolicy, class SparseMatrixPolicy>
      requires(VectorizableDense<DenseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& rate_constants,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const;

    /// @brief Returns the set of species used in a set of processes
    /// @param processes The set of processes
    /// @return The set of species used in the set of processes
    static std::set<std::string> SpeciesUsed(const std::vector<Process>& processes);
  };

  inline ProcessSet::ProcessSet(
      const std::vector<Process>& processes,
      const std::map<std::string, std::size_t>& variable_map)
      : number_of_reactants_(),
        reactant_ids_(),
        number_of_products_(),
        product_ids_(),
        yields_(),
        jacobian_process_info_(),
        jacobian_reactant_ids_(),
        jacobian_product_ids_(),
        jacobian_yields_(),
        jacobian_flat_ids_()
  {
    MICM_PROFILE_FUNCTION();

    // Set up process information for forcing calculations
    for (const auto& process : processes)
    {
      std::size_t number_of_reactants = 0;
      std::size_t number_of_products = 0;
      for (const auto& reactant : process.reactants_)
      {
        if (reactant.IsParameterized())
          continue;  // Skip reactants that are parameterizations
        if (variable_map.count(reactant.name_) < 1)
          throw std::system_error(make_error_code(MicmProcessSetErrc::ReactantDoesNotExist), reactant.name_);
        reactant_ids_.push_back(variable_map.at(reactant.name_));
        ++number_of_reactants;
      }
      for (const auto& product : process.products_)
      {
        if (product.first.IsParameterized())
          continue;  // Skip products that are parameterizations
        if (variable_map.count(product.first.name_) < 1)
          throw std::system_error(make_error_code(MicmProcessSetErrc::ProductDoesNotExist), product.first.name_);
        product_ids_.push_back(variable_map.at(product.first.name_));
        yields_.push_back(product.second);
        ++number_of_products;
      }
      number_of_reactants_.push_back(number_of_reactants);
      number_of_products_.push_back(number_of_products);
    }

    // Set up process information for Jacobian calculations
    std::vector<std::pair<std::string, std::size_t>> sorted_names(variable_map.begin(), variable_map.end());
    std::sort(sorted_names.begin(), sorted_names.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    for (const auto& independent_variable : sorted_names)
    {
      for (std::size_t i_process = 0; i_process < processes.size(); ++i_process)
      {
        const auto& process = processes[i_process];
        for (const auto& ind_reactant : process.reactants_)
        {
          if (ind_reactant.name_ != independent_variable.first)
            continue;
          ProcessInfo info;
          info.process_id_ = i_process;
          info.independent_id_ = independent_variable.second;
          info.number_of_dependent_reactants_ = 0;
          info.number_of_products_ = 0;
          bool found = false;
          for (const auto& reactant : process.reactants_)
          {
            if (reactant.IsParameterized())
              continue;  // Skip reactants that are parameterizations
            if (variable_map.count(reactant.name_) < 1)
              throw std::system_error(make_error_code(MicmProcessSetErrc::ReactantDoesNotExist), reactant.name_);
            if (variable_map.at(reactant.name_) == independent_variable.second && !found)
            {
              found = true;
              continue;
            }
            jacobian_reactant_ids_.push_back(variable_map.at(reactant.name_));
            ++info.number_of_dependent_reactants_;
          }
          for (const auto& product : process.products_)
          {
            if (product.first.IsParameterized())
              continue;  // Skip products that are parameterizations
            if (variable_map.count(product.first.name_) < 1)
              throw std::system_error(make_error_code(MicmProcessSetErrc::ProductDoesNotExist), product.first.name_);
            jacobian_product_ids_.push_back(variable_map.at(product.first.name_));
            jacobian_yields_.push_back(product.second);
            ++info.number_of_products_;
          }
          jacobian_process_info_.push_back(info);
        }
      }
    }
  };

  inline std::set<std::pair<std::size_t, std::size_t>> ProcessSet::NonZeroJacobianElements() const
  {
    MICM_PROFILE_FUNCTION();

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

  template<typename OrderingPolicy>
  inline void ProcessSet::SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix)
  {
    MICM_PROFILE_FUNCTION();

    jacobian_flat_ids_.clear();
    auto react_id = jacobian_reactant_ids_.begin();
    auto prod_id = jacobian_product_ids_.begin();
    for (const auto& process_info : jacobian_process_info_)
    {
      for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_; ++i_dep)
        jacobian_flat_ids_.push_back(matrix.VectorIndex(0, *(react_id++), process_info.independent_id_));
      jacobian_flat_ids_.push_back(matrix.VectorIndex(0, process_info.independent_id_, process_info.independent_id_));
      for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
        jacobian_flat_ids_.push_back(matrix.VectorIndex(0, *(prod_id++), process_info.independent_id_));
    }
  }

  template<typename DenseMatrixPolicy>
    requires(!VectorizableDense<DenseMatrixPolicy>)
  inline void ProcessSet::AddForcingTerms(
      const DenseMatrixPolicy& rate_constants,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    MICM_PROFILE_FUNCTION();

    // loop over grid cells
    for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
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
        {
          rate *= cell_state[react_id[i_react]];
        }

        for (std::size_t i_react = 0; i_react < number_of_reactants_[i_rxn]; ++i_react)
        {
          cell_forcing[react_id[i_react]] -= rate;
        }

        for (std::size_t i_prod = 0; i_prod < number_of_products_[i_rxn]; ++i_prod)
        {
          cell_forcing[prod_id[i_prod]] += yield[i_prod] * rate;
        }

        react_id += number_of_reactants_[i_rxn];
        prod_id += number_of_products_[i_rxn];
        yield += number_of_products_[i_rxn];
      }
    }
  };

  template<typename DenseMatrixPolicy>
    requires VectorizableDense<DenseMatrixPolicy>
  inline void ProcessSet::AddForcingTerms(
      const DenseMatrixPolicy& rate_constants,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    MICM_PROFILE_FUNCTION();

    const auto& v_rate_constants = rate_constants.AsVector();
    const auto& v_state_variables = state_variables.AsVector();
    auto& v_forcing = forcing.AsVector();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    auto v_rate_constants_begin = v_rate_constants.begin();
    // loop over all rows
    for (std::size_t i_group = 0; i_group < state_variables.NumberOfGroups(); ++i_group)
    {
      auto react_id = reactant_ids_.begin();
      auto prod_id = product_ids_.begin();
      auto yield = yields_.begin();
      const std::size_t offset_rc = i_group * rate_constants.GroupSize();
      const std::size_t offset_state = i_group * state_variables.GroupSize();
      const std::size_t offset_forcing = i_group * forcing.GroupSize();
      std::vector<double> rate(L, 0);
      const std::size_t number_of_reactions = number_of_reactants_.size();
      for (std::size_t i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
      {
        const auto v_rate_subrange_begin = v_rate_constants_begin + offset_rc + (i_rxn * L);
        rate.assign(v_rate_subrange_begin, v_rate_subrange_begin + L);
        const std::size_t number_of_reactants = number_of_reactants_[i_rxn];
        for (std::size_t i_react = 0; i_react < number_of_reactants; ++i_react)
        {
          std::size_t idx_state_variables = offset_state + react_id[i_react] * L;
          auto rate_it = rate.begin();
          auto v_state_variables_it = v_state_variables.begin() + idx_state_variables;
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(rate_it++) *= *(v_state_variables_it++);
        }
        for (std::size_t i_react = 0; i_react < number_of_reactants; ++i_react)
        {
          auto v_forcing_it = v_forcing.begin() + offset_forcing + react_id[i_react] * L;
          auto rate_it = rate.begin();
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(v_forcing_it++) -= *(rate_it++);
        }
        const std::size_t number_of_products = number_of_products_[i_rxn];
        for (std::size_t i_prod = 0; i_prod < number_of_products; ++i_prod)
        {
          auto v_forcing_it = v_forcing.begin() + offset_forcing + prod_id[i_prod] * L;
          auto rate_it = rate.begin();
          auto yield_value = yield[i_prod];
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(v_forcing_it++) += yield_value * *(rate_it++);
        }
        react_id += number_of_reactants_[i_rxn];
        prod_id += number_of_products_[i_rxn];
        yield += number_of_products_[i_rxn];
      }
    }
  }

  // Forming the Jacobian matrix "J" and returning "-J" to be consistent with the CUDA implementation
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
    requires(!VectorizableDense<DenseMatrixPolicy> || !VectorizableSparse<SparseMatrixPolicy>)
  inline void ProcessSet::SubtractJacobianTerms(
      const DenseMatrixPolicy& rate_constants,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    MICM_PROFILE_FUNCTION();

    auto cell_jacobian = jacobian.AsVector().begin();

    // loop over grid cells
    for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
    {
      auto cell_rate_constants = rate_constants[i_cell];
      auto cell_state = state_variables[i_cell];

      auto react_id = jacobian_reactant_ids_.begin();
      auto yield = jacobian_yields_.begin();
      auto flat_id = jacobian_flat_ids_.begin();

      // loop over process-dependent variable pairs
      for (const auto& process_info : jacobian_process_info_)
      {
        double d_rate_d_ind = cell_rate_constants[process_info.process_id_];
        for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          d_rate_d_ind *= cell_state[*(react_id++)];
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          cell_jacobian[*(flat_id++)] += d_rate_d_ind;
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          cell_jacobian[*(flat_id++)] -= *(yield++) * d_rate_d_ind;
      }
      // increment cell_jacobian after each grid cell
      cell_jacobian += jacobian.FlatBlockSize();
    }
  }

  // Forming the Jacobian matrix "J" and returning "-J" to be consistent with the CUDA implementation
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
    requires(VectorizableDense<DenseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  inline void ProcessSet::SubtractJacobianTerms(
      const DenseMatrixPolicy& rate_constants,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    MICM_PROFILE_FUNCTION();

    const auto& v_rate_constants = rate_constants.AsVector();
    const auto& v_state_variables = state_variables.AsVector();
    auto& v_jacobian = jacobian.AsVector();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    std::vector<double> d_rate_d_ind(L, 0);
    auto v_rate_constants_begin = v_rate_constants.begin();

    // prefetching temporary variables
    // The idea is to try to fetch the next cache line prior to operating on the current one
    // assuming the vector dimension is about the same as the cache line size
    // Adapted from https://stackoverflow.com/questions/14246976/c-how-to-force-prefetch-data-to-cache-array-loop
    double temp = 0;
    volatile double keep_temp_alive;
    constexpr std::size_t NUM_PREFETCH = 1;

    // loop over all rows
    for (std::size_t i_group = 0; i_group < state_variables.NumberOfGroups(); ++i_group)
    {
      auto react_id = jacobian_reactant_ids_.begin();
      auto yield = jacobian_yields_.begin();
      const std::size_t offset_rc = i_group * rate_constants.GroupSize();
      const std::size_t offset_state = i_group * state_variables.GroupSize();
      const std::size_t offset_jacobian = i_group * jacobian.GroupSize();
      auto flat_id = jacobian_flat_ids_.begin();

      for (const auto& process_info : jacobian_process_info_)
      {
        auto v_rate_subrange_begin = v_rate_constants_begin + offset_rc + (process_info.process_id_ * L);
        d_rate_d_ind.assign(v_rate_subrange_begin, v_rate_subrange_begin + L);
        for (std::size_t i = 0; i < std::min(NUM_PREFETCH, process_info.number_of_dependent_reactants_); ++i)
          temp += *(v_state_variables.begin() + offset_state + (react_id[i] * L));
        for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
        {
          const std::size_t idx_state_variables = offset_state + (react_id[i_react] * L);
          auto v_state_variables_it = v_state_variables.begin() + idx_state_variables;
          auto v_d_rate_d_ind_it = d_rate_d_ind.begin();
          if (i_react + NUM_PREFETCH < process_info.number_of_dependent_reactants_)
            temp += *(v_state_variables.begin() + offset_state + (react_id[i_react + NUM_PREFETCH] * L));
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(v_d_rate_d_ind_it++) *= *(v_state_variables_it++);
        }
        for (std::size_t i = 0; i < std::min(NUM_PREFETCH, process_info.number_of_dependent_reactants_); ++i)
          temp += *(v_jacobian.begin() + offset_jacobian + *(flat_id + i));
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
        {
          auto v_jacobian_it = v_jacobian.begin() + offset_jacobian + *flat_id;
          auto v_d_rate_d_ind_it = d_rate_d_ind.begin();
          if (i_dep + NUM_PREFETCH < process_info.number_of_dependent_reactants_ + 1)
            temp += *(v_jacobian.begin() + offset_jacobian + *(flat_id + NUM_PREFETCH));
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(v_jacobian_it++) += *(v_d_rate_d_ind_it++);
          ++flat_id;
        }
        for (std::size_t i = 0; i < std::min(NUM_PREFETCH, process_info.number_of_products_); ++i)
          temp += *(v_jacobian.begin() + offset_jacobian + *(flat_id + i));
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
        {
          auto v_jacobian_it = v_jacobian.begin() + offset_jacobian + *flat_id;
          auto yield_value = yield[i_dep];
          auto v_d_rate_d_ind_it = d_rate_d_ind.begin();
          if (i_dep + NUM_PREFETCH < process_info.number_of_products_)
            temp += *(v_jacobian.begin() + offset_jacobian + *(flat_id + NUM_PREFETCH));
          for (std::size_t i_cell = 0; i_cell < L; ++i_cell)
            *(v_jacobian_it++) -= yield_value * *(v_d_rate_d_ind_it++);
          ++flat_id;
        }
        react_id += process_info.number_of_dependent_reactants_;
        yield += process_info.number_of_products_;
      }
    }
    keep_temp_alive = temp;
  }

  inline std::set<std::string> ProcessSet::SpeciesUsed(const std::vector<Process>& processes)
  {
    MICM_PROFILE_FUNCTION();

    std::set<std::string> used_species;
    for (auto& process : processes)
    {
      for (auto& reactant : process.reactants_)
        used_species.insert(reactant.name_);
      for (auto& product : process.products_)
        used_species.insert(product.first.name_);
    }
    return used_species;
  }
}  // namespace micm
