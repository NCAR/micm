// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/external_model.hpp>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/solver.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <memory>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace micm
{

  /// @brief Builder of general solvers
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam RatesPolicy Calculator of forcing and Jacobian terms
  /// @tparam LinearSolverPolicy Policy for the linear solver
  /// @tparam ExternalModels Concrete external model types added via `AddExternalModel()`
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy,
      class... ExternalModels>
  class SolverBuilder
  {
    // Allow builders of any external-model pack to move state between specializations
    template<class, class, class, class, class, class, class, class...>
    friend class SolverBuilder;

   public:
    using DenseMatrixPolicyType = DenseMatrixPolicy;
    using SparseMatrixPolicyType = SparseMatrixPolicy;
    using LuDecompositionPolicyType = LuDecompositionPolicy;
    using LinearSolverPolicyType = LinearSolverPolicy;
    using StatePolicyType = StatePolicy;

   protected:
    SolverParametersPolicy options_;
    System system_;
    std::vector<Process> reactions_;
    std::vector<Constraint<DenseMatrixPolicy, SparseMatrixPolicy>> constraints_;

    std::vector<ExternalModelSystem> external_systems_;
    std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_process_sets_;
    std::vector<ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_constraints_;

    /// Owning storage of concrete external models. Moved into the built `Solver` for direct dispatch.
    std::tuple<ExternalModels...> external_models_;

    bool ignore_unused_species_ = true;
    bool reorder_state_ = true;
    bool valid_system_ = false;

   public:
    SolverBuilder() = delete;
    virtual ~SolverBuilder() = default;

    SolverBuilder(const SolverParametersPolicy& options)
      requires(sizeof...(ExternalModels) == 0)
        : options_(options)
    {
    }

    SolverBuilder(const SolverBuilder&) = default;
    SolverBuilder& operator=(const SolverBuilder&) = default;
    SolverBuilder(SolverBuilder&&) = default;
    SolverBuilder& operator=(SolverBuilder&&) = default;

   protected:
    /// Rebind constructor used by `AddExternalModel()` to seed a new specialization with
    /// state moved from a builder whose external-model pack is one element shorter.
    SolverBuilder(const SolverParametersPolicy& options, std::tuple<ExternalModels...>&& models)
        : options_(options),
          external_models_(std::move(models))
    {
    }

   public:

    /// @brief Set the chemical system
    /// @param system The chemical system
    /// @return Updated SolverBuilder
    SolverBuilder& SetSystem(const System& system)
    {
      system_ = system;
      valid_system_ = true;
      return *this;
    }

    /// @brief Set the reactions
    /// @param reactions The reactions
    /// @return Updated SolverBuilder
    SolverBuilder& SetReactions(const std::vector<Process>& reactions)
    {
      reactions_ = reactions;
      return *this;
    }

    /// @brief Set algebraic constraints for DAE solving
    /// @param constraints Vector of constraints
    /// @return Updated SolverBuilder
    SolverBuilder& SetConstraints(std::vector<Constraint<DenseMatrixPolicy, SparseMatrixPolicy>>&& constraints)
    {
      constraints_ = std::move(constraints);
      return *this;
    }

    /// @brief Set whether to ignore unused species
    /// @param ignore_unused_species True if unused species should be ignored
    /// @return Updated SolverBuilder
    SolverBuilder& SetIgnoreUnusedSpecies(bool ignore_unused_species)
    {
      ignore_unused_species_ = ignore_unused_species;
      return *this;
    }

    /// @brief Set whether to reorder the state to optimize the LU decomposition
    /// @param reorder_state True if the state should be reordered
    /// @return Updated SolverBuilder
    SolverBuilder& SetReorderState(bool reorder_state)
    {
      reorder_state_ = reorder_state;
      return *this;
    }

    /// @brief Add an external model (state variables, processes, and/or constraints)
    ///
    /// Returns a new builder whose template parameter pack is extended with `ExternalModel`.
    /// The concrete model is stored by value in a `std::tuple` that flows through `Build()`
    /// into the constructed `Solver`, which invokes the model's solve-time methods directly.
    ///
    /// If the model satisfies `HasState`, its state variables and parameters are registered.
    /// The model must satisfy at least one of `HasProcesses` or `HasConstraints`.
    template<class ExternalModel>
    auto AddExternalModel(ExternalModel model)
    {
      static_assert(
          HasProcesses<ExternalModel> || HasConstraints<ExternalModel>,
          "External model passed to AddExternalModel() must satisfy at least HasProcesses or HasConstraints");

      using NextBuilder = SolverBuilder<
          SolverParametersPolicy,
          DenseMatrixPolicy,
          SparseMatrixPolicy,
          RatesPolicy,
          LuDecompositionPolicy,
          LinearSolverPolicy,
          StatePolicy,
          ExternalModels...,
          ExternalModel>;
      auto extended_models =
          std::tuple_cat(std::move(external_models_), std::tuple<ExternalModel>(std::move(model)));
      NextBuilder next(options_, std::move(extended_models));
      next.system_ = std::move(system_);
      next.reactions_ = std::move(reactions_);
      next.constraints_ = std::move(constraints_);
      next.external_systems_ = std::move(external_systems_);
      next.external_process_sets_ = std::move(external_process_sets_);
      next.external_constraints_ = std::move(external_constraints_);
      next.ignore_unused_species_ = ignore_unused_species_;
      next.reorder_state_ = reorder_state_;
      next.valid_system_ = valid_system_;

      const auto& appended = std::get<sizeof...(ExternalModels)>(next.external_models_);
      if constexpr (HasState<ExternalModel>)
      {
        next.external_systems_.emplace_back(ExternalModelSystem{ appended });
      }
      if constexpr (HasProcesses<ExternalModel>)
      {
        next.external_process_sets_.emplace_back(
            ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ appended });
      }
      if constexpr (HasConstraints<ExternalModel>)
      {
        next.external_constraints_.emplace_back(
            ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>{ appended });
      }

      return next;
    }

    /// @brief Creates an instance of Solver with a properly configured ODE solver
    /// @return An instance of Solver
    auto Build();

   protected:
    /// @brief Returns the total state size: gas phase + all external model state variables
    Index MergedStateSize() const
    {
      Index n = system_.StateSize();
      for (const auto& m : external_systems_)
      {
        n += std::get<0>(m.state_size_func_());
      }
      return n;
    }

    /// @brief Returns all unique species names: gas phase first, then external model variables
    std::vector<std::string> MergedUniqueNames() const
    {
      auto names = system_.UniqueNames();
      for (const auto& m : external_systems_)
      {
        auto model_names = m.variable_names_func_();
        names.insert(names.end(), model_names.begin(), model_names.end());
      }
      return names;
    }

    /// @brief Checks for unused species
    /// @param rates The rates policy instance containing information about processes
    /// @throws std::system_error if an unused species is found
    void UnusedSpeciesCheck(const RatesPolicy& rates) const;

    /// @brief Gets a map of species to their index
    /// @return The species map
    std::unordered_map<std::string, Index> GetSpeciesMap() const;

    /// @brief Returns the labels of the custom parameters
    /// @return The labels of the custom parameters

    /// @brief Builds a map of unique custom parameter labels to indices by
    ///        collecting parameters from reactions, external models, and constraints.
    /// @throws MicmException if duplicate parameter labels are found.
    /// @return An unordered_map mapping each unique parameter label to its index.
    std::unordered_map<std::string, Index> GetCustomParameterMap() const;

    /// @brief Sets the absolute tolerances per species
    /// @param parameters
    /// @param species_map
    void SetAbsoluteTolerances(std::vector<Real>& tolerances, const std::unordered_map<std::string, Index>& species_map)
        const;
  };

  /// @brief Builder of CPU-based general solvers
  /// @tparam SolverParametersPolicy Parameters for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam LuDecompositionPolicy Policy for the LU decomposition
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy = Matrix<Real>,
      class SparseMatrixPolicy = SparseMatrix<Real, SparseMatrixStandardOrdering>,
      class LuDecompositionPolicy = LuDecomposition<SparseMatrixPolicy>>
  using CpuSolverBuilder = SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>,
      LuDecompositionPolicy,
      LinearSolver<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>>;

  /// @brief Builder of CPU-based general solvers with in-place LU decomposition
  /// @tparam SolverParametersPolicy Parameters for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam LuDecompositionPolicy Policy for the LU decomposition
  template<
      class SolverParametersPolicy,
      class DenseMatrix = Matrix<Real>,
      class SparseMatrixPolicy = SparseMatrix<Real, SparseMatrixStandardOrdering>,
      class LuDecompositionPolicy = LuDecompositionInPlace<SparseMatrixPolicy>>
  using CpuSolverBuilderInPlace = SolverBuilder<
      SolverParametersPolicy,
      DenseMatrix,
      SparseMatrixPolicy,
      ProcessSet<DenseMatrix, SparseMatrixPolicy>,
      LuDecompositionPolicy,
      LinearSolverInPlace<DenseMatrix, SparseMatrixPolicy, LuDecompositionPolicy>,
      State<DenseMatrix, SparseMatrixPolicy, LuDecompositionPolicy>>;

}  // namespace micm

#include "solver_builder.inl"
