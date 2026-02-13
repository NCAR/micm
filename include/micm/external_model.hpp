// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <set>

namespace micm
{
  /// @brief Represents an external model (e.g., aerosol model) that can be integrated into the MICM system
  /// @details This struct provides functions sufficient to include any state parameters and/or state variables
  ///          required by an external model to the MICM system.
  struct ExternalModelSystem
  {
    /// @brief Function to get the state size (variables, parameters) of the external model
    std::function<std::tuple<std::size_t, std::size_t>()> state_size_func_;
    /// @brief Function to get the state variable names of the external model
    std::function<std::set<std::string>()> variable_names_func_;
    /// @brief Function to get the state parameter names of the external model
    std::function<std::set<std::string>()> parameter_names_func_;


    // Default constructor
    ExternalModelSystem() = delete;
        
    /// @brief Constructor from an external model instance
    /// @tparam ModelType Type of the external model
    /// @param model Instance of the external model
    template<typename ModelType,
             typename = std::enable_if_t<!std::is_same_v<std::decay_t<ModelType>, ExternalModelSystem>>>
    ExternalModelSystem(ModelType&& model)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(std::forward<ModelType>(model));
      state_size_func_ = [shared_model]() -> std::tuple<std::size_t, std::size_t> { return shared_model->StateSize(); };
      variable_names_func_ = [shared_model]() -> std::set<std::string> { return shared_model->StateVariableNames(); };
      parameter_names_func_ = [shared_model]() -> std::set<std::string> { return shared_model->StateParameterNames(); };
    }
  };

  /// @brief Represents the processes defined by an external model (e.g., an aerosol model)
  /// @details This struct provides functions sufficient to include the processes defined by an external model
  ///          into the MICM process set, for computing forcing terms and Jacobian contributions.
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  struct ExternalModelProcessSet
  {
        /// @brief Function to get the non-zero Jacobian elements of the external model
    std::function<std::set<std::pair<std::size_t, std::size_t>>(const std::unordered_map<std::string, std::size_t>&)> non_zero_jacobian_elements_func_;
    /// @brief Function that returns the state variables that are used in the external model's processes
    std::function<std::set<std::string>()> species_used_func_;
    /// @brief Function to generate the forcing function for the external model processes
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices)>
        get_forcing_function_;
    /// @brief Function to generate the Jacobian function for the external model processes
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices,
            const SparseMatrixPolicy& jacobian)>
        get_jacobian_function_;

    /// @brief Constructor from an external model instance
    /// @tparam ModelType Type of the external model
    /// @param model Instance of the external model
    template<typename ModelType,
             typename = std::enable_if_t<!std::is_same_v<std::decay_t<ModelType>, ExternalModelProcessSet>>>
    ExternalModelProcessSet(ModelType&& model)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(std::forward<ModelType>(model));
      non_zero_jacobian_elements_func_ = [shared_model](const std::unordered_map<std::string, std::size_t>& species_map) -> std::set<std::pair<std::size_t, std::size_t>> {
          return shared_model->NonZeroJacobianElements(species_map);
      };
      species_used_func_ = [shared_model]() -> std::set<std::string> { return shared_model->SpeciesUsed(); };
      get_forcing_function_ = [shared_model](const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
                                         const std::unordered_map<std::string, std::size_t>& state_variable_indices) ->
                                         std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>
      {
        return shared_model->template ForcingFunction<DenseMatrixPolicy>(
            state_parameter_indices,
            state_variable_indices);
      };
      get_jacobian_function_ = [shared_model](const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
                const std::unordered_map<std::string, std::size_t>& state_variable_indices,
                const SparseMatrixPolicy& jacobian) ->
                std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>
      {
        return shared_model->template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
            state_parameter_indices,
            state_variable_indices,
            jacobian);
      };
    }
  };
}