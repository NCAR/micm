// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file external_model.hpp
/// @brief Defines structures to integrate external models into the MICM system
/// 
/// @section external_model_overview Overview
/// 
/// External models provide a mechanism to incorporate additional physical processes and state variables
/// into MICM that exist in separate codebases. Examples include aerosol models, cloud chemistry models,
/// precipitation, and ice models. External models can contribute both state variables (e.g., particle
/// number concentrations, condensed-phase species concentrations, droplet size distributions) and
/// processes (e.g., gas-particle partitioning, aqueous chemistry, surface reactions, evaporation/condensation).
/// 
/// @section external_model_requirements Requirements
/// 
/// To integrate an external model with MICM, the model class must implement the following interface:
/// 
/// **State Definition Methods:**
/// - `std::tuple<std::size_t, std::size_t> StateSize() const`
///   - Returns (number of state variables, number of state parameters)
/// 
/// - `std::set<std::string> StateVariableNames() const`
///   - Returns unique names for all state variables (e.g., "AEROSOL.MODE1.SO4", "CLOUD.DROPLET_NUMBER") needed
///     to describe the model's state within MICM and that will be included as ODE solver variables
/// 
/// - `std::set<std::string> StateParameterNames() const`
///   - Returns unique names for all state parameters (e.g., rate constants, physical constants) needed to
///     describe the model's state within MICM and that will be included as fixed parameters during ODE solves.
///     State parameters are almost always temperature- or pressure-dependent properties that can change between
///     solver steps but are not themselves solved for.
/// 
/// **Process Definition Methods:**
/// - `std::set<std::string> SpeciesUsed() const`
///   - Returns names of all species (gas-phase and external model variables) used in the model's processes. This
///     allows MICM to warn users if state variables are defined but not used in any processes.
/// 
/// - `std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(const std::unordered_map<std::string, std::size_t>& state_indices) const`
///   - Returns (row, column) pairs identifying non-zero elements in the Jacobian matrix, based on the model processes. The returned
///     pairs of indices are for dependent and independent variables, respectively. This allows MICM to efficiently allocate and populate the system
///     Jacobian matrix as a sparse matrix that only stores non-zero elements.
/// 
/// - `template<typename DenseMatrixPolicy> std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const`
///   - Returns a function that updates state parameters (e.g., temperature-dependent rate constants) before each solve
/// 
/// - `template<typename DenseMatrixPolicy> std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(const std::unordered_map<std::string, std::size_t>& state_parameter_indices, const std::unordered_map<std::string, std::size_t>& state_variable_indices) const`
///   - Returns a function that calculates forcing terms (tendencies) for the model's processes
///   - Function signature: `void(state_parameters, state_variables, forcing_terms)`
/// 
/// - `template<typename DenseMatrixPolicy, typename SparseMatrixPolicy> std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(const std::unordered_map<std::string, std::size_t>& state_parameter_indices, const std::unordered_map<std::string, std::size_t>& state_variable_indices, const SparseMatrixPolicy& jacobian) const`
///   - Returns a function that calculates the negative Jacobian contributions for the model's processes
///   - Function signature: `void(state_parameters, state_variables, jacobian)`
///   - Note: Each partial derivative matrix element should be subtracted from the Jacobian (i.e., `jacobian[dependent][independent] -= value`)
/// 
/// @section external_model_usage Usage
/// 
/// **1. Add the external model to the system:**
/// ```cpp
/// auto aerosol_model = MyAerosolModel(...);
/// auto system = micm::System({
///   .gas_phase_ = gas_phase,
///   .external_models_ = { aerosol_model }
/// });
/// ```
/// 
/// **2. Add the external model's processes to the solver:**
/// ```cpp
/// auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params)
///   .SetSystem(system)
///   .AddExternalModelProcesses(aerosol_model)
///   .Build();
/// ```
/// 
/// **3. Access external model state variables:**
/// ```cpp
/// auto state = solver.GetState();
/// state["AEROSOL.MODE1.SO4"] = initial_concentration;
/// double droplet_number = state["CLOUD.DROPLET_NUMBER"];
/// ```
/// 
/// @section external_model_example Example
/// 
/// See `test/integration/stub_aerosol_1.hpp` and `test/integration/stub_aerosol_2.hpp` for complete
/// example implementations demonstrating simplified gas-particle partitioning, inter-mode transfer,
/// and the use of both constant and temperature-dependent rate parameters.
/// 
#pragma once

#include <micm/system/conditions.hpp>
#include <functional>
#include <memory>
#include <string>
#include <set>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <utility>
#include <cstddef>

namespace micm
{
  /// @brief Wrapper for external model state information
  /// 
  /// This struct encapsulates an external model's state definition (variables and parameters) and provides
  /// a type-erased interface that MICM can use to query the model's state requirements. Instances are
  /// constructed automatically when an external model is added to a System via the `external_models_` field.
  /// 
  /// @note Users typically do not construct this directly; instead, pass external model instances to
  ///       `micm::System` and they will be wrapped automatically.
  struct ExternalModelSystem
  {
    /// @brief Type-erased function returning the state size (number of variables, number of parameters)
    std::function<std::tuple<std::size_t, std::size_t>()> state_size_func_;
    
    /// @brief Type-erased function returning the set of state variable names
    std::function<std::set<std::string>()> variable_names_func_;
    
    /// @brief Type-erased function returning the set of state parameter names
    std::function<std::set<std::string>()> parameter_names_func_;

    /// @brief Default constructor is deleted (must construct from an external model instance)
    ExternalModelSystem() = delete;
        
    /// @brief Constructs a type-erased wrapper from an external model instance
    /// 
    /// This constructor captures the model instance and wraps its state definition methods in
    /// std::function objects that can be called without knowledge of the original model type.
    /// 
    /// @tparam ModelType Type of the external model (must implement the external model interface)
    /// @param model External model instance (will be moved or copied into shared ownership)
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

  /// @brief Wrapper for external model process information
  /// 
  /// This struct encapsulates an external model's process definitions (forcing functions, Jacobian functions,
  /// and state parameter updates) and provides a type-erased interface that MICM can use to incorporate the
  /// model's processes into the ODE solver. Instances are constructed when `AddExternalModelProcesses()` is
  /// called on a solver builder.
  /// 
  /// The wrapped functions are used during the solve to:
  /// - Update state parameters based on environmental conditions (temperature, pressure, etc.)
  /// - Calculate forcing terms (tendencies) for the model's processes
  /// - Calculate Jacobian contributions for implicit time integration
  /// 
  /// @tparam DenseMatrixPolicy Policy for dense matrices (state variables, parameters, forcing)
  /// @tparam SparseMatrixPolicy Policy for sparse matrices (Jacobian)
  /// 
  /// @note Users typically do not construct this directly; instead, pass external model instances to
  ///       solver builder's `AddExternalModelProcesses()` method.
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  struct ExternalModelProcessSet
  {
    /// @brief Type-erased function returning non-zero Jacobian element positions
    std::function<std::set<std::pair<std::size_t, std::size_t>>(const std::unordered_map<std::string, std::size_t>&)> non_zero_jacobian_elements_func_;
    
    /// @brief Type-erased function returning the set of species used by the model's processes
    std::function<std::set<std::string>()> species_used_func_;
    
    /// @brief Type-erased function factory for state parameter updates
    /// Returns a function that updates state parameters based on environmental conditions
    std::function<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices)>
        update_state_parameters_function_;
    
    /// @brief Type-erased function factory for forcing term calculation
    /// Returns a function that computes forcing terms (tendencies) for the model's processes
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices)>
        get_forcing_function_;
    
    /// @brief Type-erased function factory for Jacobian calculation
    /// Returns a function that computes Jacobian contributions for the model's processes
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>(
            const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
            const std::unordered_map<std::string, std::size_t>& state_variable_indices,
            const SparseMatrixPolicy& jacobian)>
        get_jacobian_function_;

    /// @brief Constructs a type-erased wrapper from an external model instance
    /// 
    /// This constructor captures the model instance and wraps its process definition methods in
    /// std::function objects. The wrapped functions can then be called by MICM's solver without
    /// knowledge of the original model type, enabling seamless integration of external processes.
    /// 
    /// @tparam ModelType Type of the external model (must implement the external model interface)
    /// @param model External model instance (will be moved or copied into shared ownership)
    template<typename ModelType,
             typename = std::enable_if_t<!std::is_same_v<std::decay_t<ModelType>, ExternalModelProcessSet>>>
    ExternalModelProcessSet(ModelType&& model)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(std::forward<ModelType>(model));
      non_zero_jacobian_elements_func_ = [shared_model](const std::unordered_map<std::string, std::size_t>& species_map) -> std::set<std::pair<std::size_t, std::size_t>> {
          return shared_model->NonZeroJacobianElements(species_map);
      };
      species_used_func_ = [shared_model]() -> std::set<std::string> { return shared_model->SpeciesUsed(); };
      update_state_parameters_function_ = [shared_model](const std::unordered_map<std::string, std::size_t>& state_parameter_indices) ->
        std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>
      {
        return shared_model->template UpdateStateParametersFunction<DenseMatrixPolicy>(
            state_parameter_indices);
      };
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