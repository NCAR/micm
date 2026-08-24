// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file external_model.hpp
/// @brief Interface and setup-time wrappers for integrating external models into MICM
///
/// @section external_model_overview Overview
///
/// External models provide a mechanism to incorporate additional physical processes and state
/// variables into MICM that live in separate codebases (for treating aerosols, clouds, etc.).
/// External models can contribute state variables (e.g., condensed-phase species, droplet number),
/// state parameters (e.g., temperature-dependent rate constants), forcing/Jacobian terms, and
/// algebraic constraints (DAE mode).
///
/// @section external_model_requirements Required interface
///
/// A model integrates by exposing member functions. The solver builder detects which capability
/// groups are present with concepts and calls the corresponding methods directly at build time
/// and solve time. The concepts are:
///
/// - `HasState`         - state variables and parameters
/// - `HasProcesses`     - non-algebraic forcing and Jacobian contributions
/// - `HasConstraints`   - algebraic constraints for DAE formulation
/// - `HasInitializeConstraintParameters` - state-diagnosed constraint parameters (optional)
///
/// Every model must satisfy at least one of `HasProcesses` or `HasConstraints`. `HasState` is
/// independent and adds state variables when present.
///
/// #### State methods (satisfy `HasState`)
/// ```cpp
/// std::tuple<Index, Index> StateSize() const;
/// std::set<std::string>    StateVariableNames() const;
/// std::set<std::string>    StateParameterNames() const;
/// ```
///
/// #### Process methods (satisfy `HasProcesses`)
///
/// Build-time queries:
/// ```cpp
/// std::set<std::string>                       SpeciesUsed() const;
/// std::set<std::pair<Index, Index>>           NonZeroJacobianElements(
///     const std::unordered_map<std::string, Index>& state_indices) const;
/// ```
///
/// Build-time finalization (called once after the parameter map, species map, and Jacobian
/// sparsity are finalized so the model can cache flat indices):
/// ```cpp
/// template<class SparseMatrixPolicy>
/// void FinalizeProcessSetup(
///     const std::unordered_map<std::string, Index>& state_parameter_indices,
///     const std::unordered_map<std::string, Index>& state_variable_indices,
///     const SparseMatrixPolicy& jacobian);
/// ```
///
/// Solve-time (called by the solver each step):
/// ```cpp
/// template<class DenseMatrixPolicy>
/// void UpdateStateParameters(
///     const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
///     DenseMatrixPolicy& state_parameters) const;
///
/// template<class DenseMatrixPolicy>
/// void AddForcingTerms(
///     const DenseMatrixPolicy& state_parameters,
///     const DenseMatrixPolicy& state_variables,
///     DenseMatrixPolicy& forcing) const;
///
/// template<class DenseMatrixPolicy, class SparseMatrixPolicy>
/// void SubtractJacobianTerms(
///     const DenseMatrixPolicy& state_parameters,
///     const DenseMatrixPolicy& state_variables,
///     SparseMatrixPolicy& jacobian) const;
/// ```
///
/// Each partial derivative element should be subtracted from the Jacobian (matches the solver's
/// -J convention).
///
/// #### Constraint methods (satisfy `HasConstraints`)
///
/// Build-time queries:
/// ```cpp
/// std::set<std::string>                       ConstraintAlgebraicVariableNames() const;
/// std::set<std::string>                       ConstraintSpeciesDependencies() const;
/// std::set<std::pair<Index, Index>>           NonZeroConstraintJacobianElements(
///     const std::unordered_map<std::string, Index>& state_indices) const;
/// std::set<std::string>                       ConstraintStateParameterNames() const;
/// ```
///
/// A model may return an empty set from `ConstraintAlgebraicVariableNames()` to disable
/// constraints for the current configuration; when it does, its constraint methods are not called.
///
/// Build-time finalization:
/// ```cpp
/// template<class SparseMatrixPolicy>
/// void FinalizeConstraintSetup(
///     const std::unordered_map<std::string, Index>& state_parameter_indices,
///     const std::unordered_map<std::string, Index>& state_variable_indices,
///     const SparseMatrixPolicy& jacobian);
/// ```
///
/// Solve-time:
/// ```cpp
/// template<class DenseMatrixPolicy>
/// void UpdateConstraintStateParameters(
///     const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
///     DenseMatrixPolicy& state_parameters) const;
///
/// template<class DenseMatrixPolicy>
/// void AddConstraintResidual(
///     const DenseMatrixPolicy& state_parameters,
///     const DenseMatrixPolicy& state_variables,
///     DenseMatrixPolicy& forcing) const;
///
/// template<class DenseMatrixPolicy, class SparseMatrixPolicy>
/// void SubtractConstraintJacobian(
///     const DenseMatrixPolicy& state_parameters,
///     const DenseMatrixPolicy& state_variables,
///     SparseMatrixPolicy& jacobian) const;
/// ```
///
/// #### Initialize constraint parameters (optional, satisfy `HasInitializeConstraintParameters`)
///
/// For constraints whose parameters are diagnosed from the current state at the start of each
/// `Solve()` call (e.g., mass-conservation totals):
/// ```cpp
/// std::set<std::string> InitializeConstraintParameterNames() const;
///
/// template<class DenseMatrixPolicy>
/// void InitializeConstraintParameters(
///     const DenseMatrixPolicy& state_variables,
///     DenseMatrixPolicy& state_parameters) const;
/// ```
///
/// @section external_model_usage Usage
/// ```cpp
/// auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params)
///     .SetSystem(system)
///     .AddExternalModel(std::move(aerosol_model))
///     .Build();
///
/// auto state = solver.GetState();
/// state["AEROSOL.MODE1.SO4"] = initial_concentration;
/// ```
///
/// @section external_model_lifetime Ownership and lifetime
///
/// `SolverBuilder::AddExternalModel(m)` consumes an rvalue and returns a new builder whose
/// template parameter pack has been extended with `std::decay_t<decltype(m)>`. The concrete model
/// is stored by value in a `std::tuple` that eventually moves into the built `Solver`. The solver
/// invokes the model's solve-time methods directly (no `std::function` wrappers, no
/// `std::shared_ptr` capture), enabling inlining and safe use in `__host__ __device__` contexts.
///

#pragma once

#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <concepts>
#include <cstddef>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace micm
{
  /// @brief Type-erased build-time wrapper carrying an external model's state-definition queries
  ///
  /// Populated by `SolverBuilder::AddExternalModel()` for models that satisfy `HasState`. Only
  /// build-time queries are wrapped here; solve-time work is dispatched directly on the concrete
  /// model held in the builder / solver tuple.
  struct ExternalModelSystem
  {
    std::function<std::tuple<Index, Index>()> state_size_func_;
    std::function<std::set<std::string>()> variable_names_func_;
    std::function<std::set<std::string>()> parameter_names_func_;

    ExternalModelSystem() = delete;

    template<typename ModelType>
    ExternalModelSystem(const ModelType& model)
      requires(!std::is_same_v<std::decay_t<ModelType>, ExternalModelSystem>)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(model);
      state_size_func_ = [shared_model]() -> std::tuple<Index, Index> { return shared_model->StateSize(); };
      variable_names_func_ = [shared_model]() -> std::set<std::string> { return shared_model->StateVariableNames(); };
      parameter_names_func_ = [shared_model]() -> std::set<std::string> { return shared_model->StateParameterNames(); };
    }
  };

  /// @brief Type-erased build-time wrapper carrying an external model's process-definition queries
  ///
  /// Populated for models that satisfy `HasProcesses`. Only build-time queries (species used and
  /// Jacobian sparsity) are wrapped here. solve-time forcing/Jacobian calls are made directly on
  /// the concrete model held in the builder / solver tuple.
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  struct ExternalModelProcessSet
  {
    std::function<std::set<std::pair<Index, Index>>(const std::unordered_map<std::string, Index>&)>
        non_zero_jacobian_elements_func_;
    std::function<std::set<std::string>()> species_used_func_;

    template<typename ModelType>
    ExternalModelProcessSet(const ModelType& model)
      requires(!std::is_same_v<std::decay_t<ModelType>, ExternalModelProcessSet>)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(model);
      non_zero_jacobian_elements_func_ =
          [shared_model](const std::unordered_map<std::string, Index>& species_map) -> std::set<std::pair<Index, Index>>
      { return shared_model->NonZeroJacobianElements(species_map); };
      species_used_func_ = [shared_model]() -> std::set<std::string> { return shared_model->SpeciesUsed(); };
    }
  };

  /// @brief Type-erased build-time wrapper carrying an external model's constraint-definition queries
  ///
  /// Populated for models that satisfy `HasConstraints`. Only build-time queries are wrapped;
  /// solve-time residual/Jacobian/update calls are made directly on the concrete model held in
  /// the builder / solver tuple.
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  struct ExternalModelConstraintSet
  {
    std::function<std::set<std::string>()> algebraic_variable_names_func_;
    std::function<std::set<std::string>()> species_dependencies_func_;
    std::function<std::set<std::pair<Index, Index>>(const std::unordered_map<std::string, Index>&)>
        non_zero_jacobian_elements_func_;
    std::function<std::set<std::string>()> state_parameter_names_func_;
    /// Empty when the model does not need state-diagnosed constraint parameters
    std::function<std::set<std::string>()> initialize_constraint_parameter_names_func_;

    template<typename ModelType>
    ExternalModelConstraintSet(const ModelType& model)
      requires(!std::is_same_v<std::decay_t<ModelType>, ExternalModelConstraintSet>)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(model);
      algebraic_variable_names_func_ = [shared_model]() -> std::set<std::string>
      { return shared_model->ConstraintAlgebraicVariableNames(); };
      species_dependencies_func_ = [shared_model]() -> std::set<std::string>
      { return shared_model->ConstraintSpeciesDependencies(); };
      non_zero_jacobian_elements_func_ =
          [shared_model](const std::unordered_map<std::string, Index>& species_map) -> std::set<std::pair<Index, Index>>
      { return shared_model->NonZeroConstraintJacobianElements(species_map); };
      state_parameter_names_func_ = [shared_model]() -> std::set<std::string>
      { return shared_model->ConstraintStateParameterNames(); };

      if constexpr (requires { shared_model->InitializeConstraintParameterNames(); })
      {
        initialize_constraint_parameter_names_func_ = [shared_model]() -> std::set<std::string>
        { return shared_model->InitializeConstraintParameterNames(); };
      }
      else
      {
        initialize_constraint_parameter_names_func_ = []() -> std::set<std::string> { return {}; };
      }
    }
  };

  /// @brief Model provides state variables and parameters
  template<typename T>
  concept HasState = requires(const T& m) {
    { m.StateSize() } -> std::same_as<std::tuple<Index, Index>>;
    { m.StateVariableNames() } -> std::same_as<std::set<std::string>>;
    { m.StateParameterNames() } -> std::same_as<std::set<std::string>>;
  };

  /// @brief Model provides forcing/Jacobian process contributions
  template<typename T>
  concept HasProcesses = requires(const T& m, const std::unordered_map<std::string, Index>& map) {
    { m.SpeciesUsed() } -> std::same_as<std::set<std::string>>;
    { m.NonZeroJacobianElements(map) } -> std::same_as<std::set<std::pair<Index, Index>>>;
  };

  /// @brief Model provides algebraic constraint contributions
  template<typename T>
  concept HasConstraints = requires(const T& m) {
    { m.ConstraintAlgebraicVariableNames() } -> std::same_as<std::set<std::string>>;
    { m.ConstraintSpeciesDependencies() } -> std::same_as<std::set<std::string>>;
  };

  /// @brief Model provides state-diagnosed constraint parameters (called once per Solve())
  template<typename T>
  concept HasInitializeConstraintParameters = requires(const T& m) {
    { m.InitializeConstraintParameterNames() } -> std::same_as<std::set<std::string>>;
  };
}  // namespace micm
