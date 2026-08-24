// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/external_model.hpp>
#include <micm/process/process.hpp>
#include <micm/util/error.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace micm
{

  /// @brief Solver function calculators for a collection of processes
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  class ProcessSet
  {
    using DenseMatrix = DenseMatrixPolicy;
    using SparseMatrix = SparseMatrixPolicy;
    template<class U>
    using Scalar = typename SparseMatrix::template ScalarType<U>;
    template<class U>
    using Vector = typename SparseMatrix::template VectorType<U>;
    template<class U>
    using VectorView = typename SparseMatrix::template VectorType<U>::ConstViewType;

   public:
    /// @brief Process information for use in setting Jacobian elements
    struct ProcessInfo
    {
      Index process_id_;
      Index independent_id_;
      Index number_of_dependent_reactants_;
      Index number_of_products_;
    };
    struct Views
    {
      VectorView<Index> number_of_reactants_;
      VectorView<Index> reactant_ids_;
      VectorView<Index> number_of_products_;
      VectorView<Index> product_ids_;
      VectorView<Real> yields_;
      VectorView<ProcessInfo> jacobian_process_info_;
      VectorView<Index> jacobian_reactant_ids_;
      VectorView<Index> jacobian_product_ids_;
      VectorView<Real> jacobian_yields_;
      VectorView<Index> jacobian_flat_ids_;
      VectorView<Bool> is_algebraic_variable_;

      Views() = default;

      Views(
          const Vector<Index>& number_of_reactants,
          const Vector<Index>& reactant_ids,
          const Vector<Index>& number_of_products,
          const Vector<Index>& product_ids,
          const Vector<Real>& yields,
          const Vector<ProcessInfo>& jacobian_process_info,
          const Vector<Index>& jacobian_reactant_ids,
          const Vector<Index>& jacobian_product_ids,
          const Vector<Real>& jacobian_yields,
          const Vector<Index>& jacobian_flat_ids,
          const Vector<Bool>& is_algebraic_variable)
          : number_of_reactants_(number_of_reactants.GetView()),
            reactant_ids_(reactant_ids.GetView()),
            number_of_products_(number_of_products.GetView()),
            product_ids_(product_ids.GetView()),
            yields_(yields.GetView()),
            jacobian_process_info_(jacobian_process_info.GetView()),
            jacobian_reactant_ids_(jacobian_reactant_ids.GetView()),
            jacobian_product_ids_(jacobian_product_ids.GetView()),
            jacobian_yields_(jacobian_yields.GetView()),
            jacobian_flat_ids_(jacobian_flat_ids.GetView()),
            is_algebraic_variable_(is_algebraic_variable.GetView())
      {
      }
    };

   protected:
    Vector<Index> number_of_reactants_;
    Vector<Index> reactant_ids_;
    Vector<Index> number_of_products_;
    Vector<Index> product_ids_;
    Vector<Real> yields_;
    Vector<ProcessInfo> jacobian_process_info_;
    Vector<Index> jacobian_reactant_ids_;
    Vector<Index> jacobian_product_ids_;
    Vector<Real> jacobian_yields_;
    Vector<Index> jacobian_flat_ids_;
    Vector<Bool> is_algebraic_variable_;
    Views views_;

    std::unordered_map<std::string, Index> variable_map_;

    std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_process_sets_;
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>>
        external_forcing_functions_;
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>>
        external_jacobian_functions_;

   public:
    /// @brief Default constructor
    ProcessSet() = default;

    ProcessSet(ProcessSet&& other) noexcept;
    ProcessSet& operator=(ProcessSet&& other) noexcept;

    /// @brief Constructs a ProcessSet by mapping species in each process to their corresponding indices
    ///        Initializes internal data structures related to a set of processes, mapping them to variable indices
    ///        using a provided variable_map. Also prepares the data needed for computing Jacobian contributions.
    /// @param processes A list of processes, each with reactants and products
    /// @param variable_map A map from species names to their corresponding index in the solver's state
    /// @throws std::system_error If a reactant or product name in a process is not found in variable_map
    ProcessSet(const std::vector<Process>& processes, const std::unordered_map<std::string, Index>& variable_map);

    /// @brief Constructs a ProcessSet as above, but also includes contributions from external models
    /// @param processes A list of processes, each with reactants and products
    /// @param variable_map A map from species names to their corresponding index in the solver's state
    /// @param external_process_sets A list of external process sets that provide additional processes and Jacobian
    /// contributions
    /// @throws std::system_error If a reactant or product name in a process is not found in variable_map
    ProcessSet(
        const std::vector<Process>& processes,
        const std::unordered_map<std::string, Index>& variable_map,
        const std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>>& external_process_sets)
        : ProcessSet(processes, variable_map)
    {
      external_process_sets_ = external_process_sets;
    }

    virtual ~ProcessSet() = default;

    /// @brief Returns the positions of all non-zero Jacobian elements
    /// @return A set of (row, column) index pairs, each representing a non-zero entry
    std::set<std::pair<Index, Index>> NonZeroJacobianElements() const;

    /// @brief Computes and stores flat (1D) indices for non-zero Jacobian elements
    ///        Stores combination of process ids and reactant ids to support column-wise Jacobian updates.
    /// @param matrix The sparse Jacobian matrix used to compute flat indices.
    void SetJacobianFlatIds(const SparseMatrixPolicy& matrix);

    /// @brief Marks species rows that should be treated as algebraic (constraints replace ODE rows)
    /// @param variable_ids Set of variable ids whose forcing/Jacobian rows should not receive kinetic contributions
    void SetAlgebraicVariableIds(const std::set<Index>& variable_ids);

    /// @brief Sets external model functions for forcing terms and Jacobian contributions
    /// @param state_parameter_indices Map of state parameter names to their indices
    /// @param state_variable_indices Map of state variable names to their indices
    /// @param jacobian The sparse Jacobian matrix used by the solver
    void SetExternalModelFunctions(
        const std::unordered_map<std::string, Index>& state_parameter_indices,
        const std::unordered_map<std::string, Index>& state_variable_indices,
        const SparseMatrixPolicy& jacobian);

    /// @brief Adds forcing terms for the set of processes for the current conditions
    /// @param state Current state containing rate constants and other relevant data
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param forcing Forcing terms for each state variable (grid cell, state variable)
    template<class StatePolicy>
    void AddForcingTerms(const StatePolicy& state, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing)
        const;

    /// @brief Subtracts Jacobian terms for the set of processes for the current conditions
    /// @param state Current state containing rate constants and other relevant data
    /// @param state_variables Current state variable values (grid cell, state variable)
    /// @param jacobian Jacobian matrix for the system (grid cell, dependent variable, independent variable)
    template<class StatePolicy>
    void SubtractJacobianTerms(
        const StatePolicy& state,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const;

    /// @brief Extracts all species involved in the given processes
    /// @param processes A list of Process objects, each with reactants and products
    /// @return A set of species names
    std::set<std::string> SpeciesUsed(const std::vector<Process>& processes) const;
  };

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::ProcessSet(ProcessSet&& other) noexcept
      : number_of_reactants_(std::move(other.number_of_reactants_)),
        reactant_ids_(std::move(other.reactant_ids_)),
        number_of_products_(std::move(other.number_of_products_)),
        product_ids_(std::move(other.product_ids_)),
        yields_(std::move(other.yields_)),
        jacobian_process_info_(std::move(other.jacobian_process_info_)),
        jacobian_reactant_ids_(std::move(other.jacobian_reactant_ids_)),
        jacobian_product_ids_(std::move(other.jacobian_product_ids_)),
        jacobian_yields_(std::move(other.jacobian_yields_)),
        jacobian_flat_ids_(std::move(other.jacobian_flat_ids_)),
        is_algebraic_variable_(std::move(other.is_algebraic_variable_)),
        views_(
            number_of_reactants_,
            reactant_ids_,
            number_of_products_,
            product_ids_,
            yields_,
            jacobian_process_info_,
            jacobian_reactant_ids_,
            jacobian_product_ids_,
            jacobian_yields_,
            jacobian_flat_ids_,
            is_algebraic_variable_),
        variable_map_(std::move(other.variable_map_)),
        external_process_sets_(std::move(other.external_process_sets_)),
        external_forcing_functions_(std::move(other.external_forcing_functions_)),
        external_jacobian_functions_(std::move(other.external_jacobian_functions_))
  {
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>& ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::operator=(
      ProcessSet&& other) noexcept
  {
    if (this != &other)
    {
      number_of_reactants_ = std::move(other.number_of_reactants_);
      reactant_ids_ = std::move(other.reactant_ids_);
      number_of_products_ = std::move(other.number_of_products_);
      product_ids_ = std::move(other.product_ids_);
      yields_ = std::move(other.yields_);
      jacobian_process_info_ = std::move(other.jacobian_process_info_);
      jacobian_reactant_ids_ = std::move(other.jacobian_reactant_ids_);
      jacobian_product_ids_ = std::move(other.jacobian_product_ids_);
      jacobian_yields_ = std::move(other.jacobian_yields_);
      jacobian_flat_ids_ = std::move(other.jacobian_flat_ids_);
      is_algebraic_variable_ = std::move(other.is_algebraic_variable_);
      views_ = Views(
          number_of_reactants_,
          reactant_ids_,
          number_of_products_,
          product_ids_,
          yields_,
          jacobian_process_info_,
          jacobian_reactant_ids_,
          jacobian_product_ids_,
          jacobian_yields_,
          jacobian_flat_ids_,
          is_algebraic_variable_);
      variable_map_ = std::move(other.variable_map_);
      external_process_sets_ = std::move(other.external_process_sets_);
      external_forcing_functions_ = std::move(other.external_forcing_functions_);
      external_jacobian_functions_ = std::move(other.external_jacobian_functions_);
    }
    return *this;
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::ProcessSet(
      const std::vector<Process>& processes,
      const std::unordered_map<std::string, Index>& variable_map)
      : number_of_reactants_(),
        reactant_ids_(),
        number_of_products_(),
        product_ids_(),
        yields_(),
        jacobian_process_info_(),
        jacobian_reactant_ids_(),
        jacobian_product_ids_(),
        jacobian_yields_(),
        jacobian_flat_ids_(),
        is_algebraic_variable_(variable_map.size(), micm::Bool(false)),
        views_(),
        variable_map_(variable_map)
  {
    // For each process, look up each reactant name in variable_map and
    // store the corresponding index
    std::vector<Index> number_of_reactants_temp;
    std::vector<Index> reactant_ids_temp;
    std::vector<Index> number_of_products_temp;
    std::vector<Index> product_ids_temp;
    std::vector<Real> yields_temp;
    for (const auto& process : processes)
    {
      const auto& reaction = process.process_;
      Index number_of_reactants = 0;
      Index number_of_products = 0;
      for (const auto& reactant : reaction.reactants_)
      {
        if (reactant.IsParameterized())
        {
          continue;  // Skip reactants that are parameterizations
        }
        if (!variable_map.contains(reactant.name_))
        {
          throw MicmException(MICM_ERROR_CATEGORY_PROCESS, MICM_PROCESS_ERROR_CODE_REACTANT_DOES_NOT_EXIST, reactant.name_);
        }
        reactant_ids_temp.push_back(variable_map.at(reactant.name_));
        ++number_of_reactants;
      }
      // Store product indices and yields
      for (const auto& product : reaction.products_)
      {
        if (product.species_.IsParameterized())
        {
          continue;  // Skip products that are parameterizations
        }
        if (!variable_map.contains(product.species_.name_))
        {
          throw MicmException(
              MICM_ERROR_CATEGORY_PROCESS, MICM_PROCESS_ERROR_CODE_PRODUCT_DOES_NOT_EXIST, product.species_.name_);
        }
        product_ids_temp.push_back(variable_map.at(product.species_.name_));
        yields_temp.push_back(product.coefficient_);
        ++number_of_products;
      }
      // Record how many reactants and products were processed for each process
      number_of_reactants_temp.push_back(number_of_reactants);
      number_of_products_temp.push_back(number_of_products);
    }

    // Set up process information for Jacobian calculations.
    //
    // For every independent variable (species) used as a reactant in a process, create a
    // ProcessInfo record. Rather than scanning every process for every species
    // (O(species x reactions), which is prohibitive for large mechanisms), collect
    // (independent_id, process_id) jobs in a SINGLE pass over the processes and then
    // stable-sort them by independent_id. The emitted records are then ordered by
    // independent variable and then by process which is identical to the species-by-species
    // scan, but the cost is O(reactions x reactants + J log J) instead of
    // O(species x reactions).
    std::vector<std::pair<Index, Index>> jacobian_columns;  // (independent_id, i_process)
    std::vector<ProcessInfo> jacobian_process_info_temp;
    std::vector<Index> jacobian_reactant_ids_temp;
    std::vector<Index> jacobian_product_ids_temp;
    std::vector<Real> jacobian_yields_temp;
    for (Index i_process = 0; i_process < processes.size(); ++i_process)
    {
      const auto& reaction = processes[i_process].process_;
      for (const auto& ind_reactant : reaction.reactants_)
      {
        auto it = variable_map.find(ind_reactant.name_);
        if (it == variable_map.end())
        {
          continue;
        }
        jacobian_columns.emplace_back(it->second, i_process);
      }
    }
    std::stable_sort(
        jacobian_columns.begin(), jacobian_columns.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    for (const auto& column : jacobian_columns)
    {
      const Index independent_id = column.first;
      const Index i_process = column.second;
      const auto& reaction = processes[i_process].process_;
      ProcessInfo info;
      info.process_id_ = i_process;
      info.independent_id_ = independent_id;
      info.number_of_dependent_reactants_ = 0;
      info.number_of_products_ = 0;

      // Collect other (dependent) reactants and products
      // because our reactants can be duplicated, (2B could be 2B or B + B), we need to detect when we've already seen a
      // reactant
      bool found = false;
      for (const auto& reactant : reaction.reactants_)
      {
        if (reactant.IsParameterized())
        {
          continue;  // Skip reactants that are parameterizations
        }
        const Index id = variable_map.at(reactant.name_);
        if (id == independent_id && !found)
        {
          found = true;
          continue;
        }
        jacobian_reactant_ids_temp.push_back(id);
        ++info.number_of_dependent_reactants_;
      }
      for (const auto& product : reaction.products_)
      {
        if (product.species_.IsParameterized())
        {
          continue;  // Skip products that are parameterizations
        }
        jacobian_product_ids_temp.push_back(variable_map.at(product.species_.name_));
        jacobian_yields_temp.push_back(product.coefficient_);
        ++info.number_of_products_;
      }
      jacobian_process_info_temp.push_back(info);
    }
    number_of_reactants_ = number_of_reactants_temp;
    reactant_ids_ = reactant_ids_temp;
    number_of_products_ = number_of_products_temp;
    product_ids_ = product_ids_temp;
    yields_ = yields_temp;
    jacobian_process_info_ = jacobian_process_info_temp;
    jacobian_reactant_ids_ = jacobian_reactant_ids_temp;
    jacobian_product_ids_ = jacobian_product_ids_temp;
    jacobian_yields_ = jacobian_yields_temp;
    number_of_reactants_.CopyToDevice();
    reactant_ids_.CopyToDevice();
    number_of_products_.CopyToDevice();
    product_ids_.CopyToDevice();
    yields_.CopyToDevice();
    jacobian_process_info_.CopyToDevice();
    jacobian_reactant_ids_.CopyToDevice();
    jacobian_product_ids_.CopyToDevice();
    jacobian_yields_.CopyToDevice();
    views_ = Views(
        number_of_reactants_,
        reactant_ids_,
        number_of_products_,
        product_ids_,
        yields_,
        jacobian_process_info_,
        jacobian_reactant_ids_,
        jacobian_product_ids_,
        jacobian_yields_,
        jacobian_flat_ids_,
        is_algebraic_variable_);
  };

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline std::set<std::pair<Index, Index>> ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::NonZeroJacobianElements() const
  {
    std::set<std::pair<Index, Index>> ids;
    auto react_id = reactant_ids_.begin();
    auto prod_id = product_ids_.begin();
    for (Index i_rxn = 0; i_rxn < number_of_reactants_.size(); ++i_rxn)
    {
      for (Index i_ind = 0; i_ind < number_of_reactants_[i_rxn]; ++i_ind)
      {
        // For each reactant, collect the Jacobian contributing indices
        for (Index i_dep = 0; i_dep < number_of_reactants_[i_rxn]; ++i_dep)
        {
          ids.insert(std::make_pair(react_id[i_dep], react_id[i_ind]));
        }
        // For each product, collect the Jacobian contributing indices
        for (Index i_dep = 0; i_dep < number_of_products_[i_rxn]; ++i_dep)
        {
          ids.insert(std::make_pair(prod_id[i_dep], react_id[i_ind]));
        }
      }
      // Adavance iterators using the number of reactants/products in each process
      react_id += number_of_reactants_[i_rxn];
      prod_id += number_of_products_[i_rxn];
    }

    // Add Jacobian elements from external process sets
    for (const auto& process_set : external_process_sets_)
    {
      auto external_jac_elements = process_set.non_zero_jacobian_elements_func_(variable_map_);
      ids.insert(external_jac_elements.begin(), external_jac_elements.end());
    }
    return ids;
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline void ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetJacobianFlatIds(const SparseMatrixPolicy& matrix)
  {
    std::vector<Index> jacobian_flat_ids_temp;
    jacobian_flat_ids_temp.clear();
    jacobian_flat_ids_temp.reserve(
        jacobian_reactant_ids_.size() + jacobian_process_info_.size() + jacobian_product_ids_.size());
    auto react_id = jacobian_reactant_ids_.begin();
    auto prod_id = jacobian_product_ids_.begin();
    // Algebraic rows may be pruned from sparsity; keep placeholder ids so the update loops stay aligned.
    constexpr Index skipped_flat_id = 0;
    for (const auto& process_info : jacobian_process_info_)
    {
      for (Index i_dep = 0; i_dep < process_info.number_of_dependent_reactants_; ++i_dep)
      {
        const Index row_id = *(react_id++);
        jacobian_flat_ids_temp.push_back(
            is_algebraic_variable_[row_id] ? skipped_flat_id : matrix.VectorIndex(0, row_id, process_info.independent_id_));
      }
      jacobian_flat_ids_temp.push_back(
          is_algebraic_variable_[process_info.independent_id_]
              ? skipped_flat_id
              : matrix.VectorIndex(0, process_info.independent_id_, process_info.independent_id_));
      for (Index i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
      {
        const Index row_id = *(prod_id++);
        jacobian_flat_ids_temp.push_back(
            is_algebraic_variable_[row_id] ? skipped_flat_id : matrix.VectorIndex(0, row_id, process_info.independent_id_));
      }
    }
    jacobian_flat_ids_ = jacobian_flat_ids_temp;
    jacobian_flat_ids_.CopyToDevice();
    views_.jacobian_flat_ids_ = std::as_const(jacobian_flat_ids_).GetView();
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline void ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetAlgebraicVariableIds(const std::set<Index>& variable_ids)
  {
    std::fill(is_algebraic_variable_.begin(), is_algebraic_variable_.end(), false);
    for (const auto variable_id : variable_ids)
    {
      is_algebraic_variable_[variable_id] = micm::Bool(true);
    }
    is_algebraic_variable_.CopyToDevice();
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline void ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetExternalModelFunctions(
      const std::unordered_map<std::string, Index>& state_parameter_indices,
      const std::unordered_map<std::string, Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    external_forcing_functions_.clear();
    external_jacobian_functions_.clear();
    for (const auto& process_set : external_process_sets_)
    {
      external_forcing_functions_.push_back(
          process_set.get_forcing_function_(state_parameter_indices, state_variable_indices));
      external_jacobian_functions_.push_back(
          process_set.get_jacobian_function_(state_parameter_indices, state_variable_indices, jacobian));
    }
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  template<class StatePolicy>
  inline void ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::AddForcingTerms(
      const StatePolicy& state,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const Index n_rxn = number_of_reactants_.size();
    const auto& views = views_;
    const DenseMatrix& rate_constants = state.rate_constants_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrix::ViewType& forcing_view,
            const typename DenseMatrix::ConstViewType& state_view,
            const typename DenseMatrix::ConstViewType& rc_view) {
          auto react_id = views.reactant_ids_.begin();
          auto prod_id = views.product_ids_.begin();
          auto yield = views.yields_.begin();
          auto rate = forcing_view.GetRowVariable();
          for (Index i_rxn = 0; i_rxn < n_rxn; ++i_rxn)
          {
            // Calculate the rate as the rate constant times the product of all reactant concentrations
            rc_view.Copy(rate, rc_view.GetConstColumnView(i_rxn));
            for (Index i_react = 0; i_react < views.number_of_reactants_[i_rxn]; ++i_react)
            {
              state_view.ForEachRow(
                  [](Real& rate, const Real& reactant) { rate *= reactant; },
                  rate,
                  state_view.GetConstColumnView(react_id[i_react]));
            }
            // Subtract the rate from reactant forcings
            for (Index i_react = 0; i_react < views.number_of_reactants_[i_rxn]; ++i_react)
            {
              const Index row_id = react_id[i_react];
              if (!views.is_algebraic_variable_[row_id])
              {
                forcing_view.ForEachRow(
                    [](Real& forcing, const Real& rate) { forcing -= rate; }, forcing_view.GetColumnView(row_id), rate);
              }
            }
            // Add the rate (scaled by yield) to the product forcings
            for (Index i_prod = 0; i_prod < views.number_of_products_[i_rxn]; ++i_prod)
            {
              const Index row_id = prod_id[i_prod];
              if (!views.is_algebraic_variable_[row_id])
              {
                const Real prod_yield = yield[i_prod];
                forcing_view.ForEachRow(
                    [&prod_yield](Real& forcing, const Real& rate) { forcing += prod_yield * rate; },
                    forcing_view.GetColumnView(row_id),
                    rate);
              }
            }
            // Update iterators based on how many reactants/products each reaction has
            react_id += views.number_of_reactants_[i_rxn];
            prod_id += views.number_of_products_[i_rxn];
            yield += views.number_of_products_[i_rxn];
          }
        },
        forcing,
        state_variables,
        rate_constants)(forcing, state_variables, rate_constants);

    // Add forcing contributions from external models
    for (const auto& add_forcing_function : external_forcing_functions_)
    {
      add_forcing_function(state.custom_rate_parameters_, state_variables, forcing);
    }
  };

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  template<class StatePolicy>
  inline void ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SubtractJacobianTerms(
      const StatePolicy& state,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    const Index n_rxn = number_of_reactants_.size();
    const auto& views = views_;
    const DenseMatrix& rate_constants = state.rate_constants_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename SparseMatrix::ViewType& jacobian_view,
            const typename DenseMatrix::ConstViewType& state_view,
            const typename DenseMatrix::ConstViewType& rc_view) {
          auto react_id = views.jacobian_reactant_ids_.begin();
          auto prod_id = views.jacobian_product_ids_.begin();
          auto yield = views.jacobian_yields_.begin();
          auto flat_id = views.jacobian_flat_ids_.begin();
          auto d_rate_d_ind = jacobian_view.GetBlockVariable();
          // loop over process-dependent variable pairs
          for (const auto& process_info : views.jacobian_process_info_)
          {
            rc_view.Copy(d_rate_d_ind, rc_view.GetConstColumnView(process_info.process_id_));
            for (Index i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
            {
              state_view.ForEachRow(
                  [](Real& dr_di, const Real& reactant) { dr_di *= reactant; },
                  d_rate_d_ind,
                  state_view.GetConstColumnView(react_id[i_react]));
            }
            for (Index i_dep = 0; i_dep < process_info.number_of_dependent_reactants_; ++i_dep)
            {
              const Index row_id = react_id[i_dep];
              if (!views.is_algebraic_variable_[row_id])
              {
                jacobian_view.ForEachBlock(
                    [](Real& jacobian, const Real& dr_di) { jacobian += dr_di; },
                    jacobian_view.GetBlockView(*flat_id),
                    d_rate_d_ind);
              }
              ++flat_id;
            }
            if (!views.is_algebraic_variable_[process_info.independent_id_])
            {
              jacobian_view.ForEachBlock(
                  [](Real& jacobian, const Real& dr_di) { jacobian += dr_di; },
                  jacobian_view.GetBlockView(*flat_id),
                  d_rate_d_ind);
            }
            ++flat_id;
            for (Index i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
            {
              const Index row_id = prod_id[i_dep];
              if (!views.is_algebraic_variable_[row_id])
              {
                const Real prod_yield = yield[i_dep];
                jacobian_view.ForEachBlock(
                    [prod_yield](Real& jacobian, const Real& dr_di) { jacobian -= prod_yield * dr_di; },
                    jacobian_view.GetBlockView(*flat_id),
                    d_rate_d_ind);
              }
              ++flat_id;
            }
            react_id += process_info.number_of_dependent_reactants_;
            prod_id += process_info.number_of_products_;
            yield += process_info.number_of_products_;
          }
        },
        jacobian,
        state_variables,
        rate_constants)(jacobian, state_variables, rate_constants);

    // Add Jacobian contributions from external models
    for (const auto& add_jacobian_function : external_jacobian_functions_)
    {
      add_jacobian_function(state.custom_rate_parameters_, state_variables, jacobian);
    }
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  inline std::set<std::string> ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SpeciesUsed(
      const std::vector<Process>& processes) const
  {
    std::set<std::string> used_species;
    for (const auto& process : processes)
    {
      const auto& reaction = process.process_;
      for (const auto& reactant : reaction.reactants_)
      {
        used_species.insert(reactant.name_);
      }
      for (const auto& product : reaction.products_)
      {
        used_species.insert(product.species_.name_);
      }
    }

    // Include species used in external process sets
    for (const auto& process_set : external_process_sets_)
    {
      auto external_species_used = process_set.species_used_func_();
      used_species.insert(external_species_used.begin(), external_species_used.end());
    }

    return used_species;
  }
}  // namespace micm
