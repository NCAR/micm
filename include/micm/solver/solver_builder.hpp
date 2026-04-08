// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
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
#include <micm/util/vector_matrix.hpp>

#include <memory>
#include <sstream>
#include <unordered_map>

namespace micm
{

  /// @brief Builder of general solvers
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam RatesPolicy Calculator of forcing and Jacobian terms
  /// @tparam LinearSolverPolicy Policy for the linear solver
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  class SolverBuilder
  {
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
    std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_models_;
    std::vector<ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_constraint_models_;
    std::vector<Constraint> constraints_;
    bool ignore_unused_species_ = true;
    bool reorder_state_ = true;
    bool valid_system_ = false;

   public:
    SolverBuilder() = delete;
    virtual ~SolverBuilder() = default;

    SolverBuilder(const SolverParametersPolicy& options)
        : options_(options)
    {
    }

    // Copy constructor
    SolverBuilder(const SolverBuilder& other) = default;

    // Copy assignment
    SolverBuilder& operator=(const SolverBuilder& other) = default;

    // Default move operations
    SolverBuilder(SolverBuilder&&) = default;
    SolverBuilder& operator=(SolverBuilder&&) = default;

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

    /// @brief Add processes from an external model
    /// @param model The external model
    /// @return Updated SolverBuilder
    /// @deprecated Use AddExternalModel() instead
    template<class ExternalModel>
    SolverBuilder& AddExternalModelProcesses(ExternalModel&& model)
    {
      external_models_.emplace_back(
          ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::forward<decltype(model)>(model) });
      return *this;
    }

    /// @brief Add constraints from an external model
    ///
    /// Only wraps constraint information. The model must satisfy the HasConstraints concept.
    /// Use this when processes are added separately or not needed.
    /// @param model The external model
    /// @return Updated SolverBuilder
    template<class ExternalModel>
    SolverBuilder& AddExternalModelConstraints(ExternalModel&& model)
    {
      static_assert(
          HasConstraints<std::decay_t<ExternalModel>>,
          "External model passed to AddExternalModelConstraints() must satisfy the HasConstraints concept");
      external_constraint_models_.emplace_back(
          ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::forward<ExternalModel>(model) });
      return *this;
    }

    /// @brief Add an external model (processes and optionally constraints)
    ///
    /// If the model satisfies the HasConstraints concept, both process and constraint
    /// wrappers are created. Otherwise, only processes are wrapped.
    /// @param model The external model
    /// @return Updated SolverBuilder
    template<class ExternalModel>
    SolverBuilder& AddExternalModel(ExternalModel&& model)
    {
      // Always wrap process info
      auto model_copy = model;
      external_models_.emplace_back(
          ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::move(model_copy) });
      // Conditionally wrap constraint info
      if constexpr (HasConstraints<std::decay_t<ExternalModel>>)
      {
        external_constraint_models_.emplace_back(
            ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::forward<ExternalModel>(model) });
      }
      return *this;
    }

    /// @brief Set algebraic constraints for DAE solving
    /// @param constraints Vector of constraints
    /// @return Updated SolverBuilder
    SolverBuilder& SetConstraints(std::vector<Constraint>&& constraints)
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

    /// @brief Creates an instance of Solver with a properly configured ODE solver
    /// @return An instance of Solver
    auto Build();

   protected:
    /// @brief Checks for unused species
    /// @param rates The rates policy instance containing information about processes
    /// @throws std::system_error if an unused species is found
    void UnusedSpeciesCheck(const RatesPolicy& rates) const;

    /// @brief Gets a map of species to their index
    /// @return The species map
    std::unordered_map<std::string, std::size_t> GetSpeciesMap() const;

    /// @brief Returns the labels of the custom parameters
    /// @return The labels of the custom parameters

    /// @brief Builds a map of unique custom parameter labels to indices by
    ///        collecting parameters from reactions, external models, and constraints.
    /// @throws MicmException if duplicate parameter labels are found.
    /// @return An unordered_map mapping each unique parameter label to its index.
    std::unordered_map<std::string, std::size_t> GetCustomParameterMap() const;

    /// @brief Sets the absolute tolerances per species
    /// @param parameters
    /// @param species_map
    void SetAbsoluteTolerances(
        std::vector<double>& tolerances,
        const std::unordered_map<std::string, std::size_t>& species_map) const;
  };

  /// @brief Builder of CPU-based general solvers
  /// @tparam SolverParametersPolicy Parameters for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam LuDecompositionPolicy Policy for the LU decomposition
  /// @tparam LMatrixPolicy Policy for the Lower matrix
  /// @tparam UMatrixPolicy Policy for the Upper matrix
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy = Matrix<double>,
      class SparseMatrixPolicy = SparseMatrix<double, SparseMatrixStandardOrdering>,
      class LuDecompositionPolicy = LuDecomposition,
      class LMatrixPolicy = SparseMatrixPolicy,
      class UMatrixPolicy = SparseMatrixPolicy>
  using CpuSolverBuilder = SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>,
      LuDecompositionPolicy,
      LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>>;

  /// @brief Builder of CPU-based general solvers with in-place LU decomposition
  /// @tparam SolverParametersPolicy Parameters for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  /// @tparam LuDecompositionPolicy Policy for the LU decomposition
  template<
      class SolverParametersPolicy,
      class DenseMatrix = Matrix<double>,
      class SparseMatrixPolicy = SparseMatrix<double, SparseMatrixStandardOrdering>,
      class LuDecompositionPolicy = LuDecompositionInPlace>
  using CpuSolverBuilderInPlace = SolverBuilder<
      SolverParametersPolicy,
      DenseMatrix,
      SparseMatrixPolicy,
      ProcessSet<DenseMatrix, SparseMatrixPolicy>,
      LuDecompositionPolicy,
      LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>,
      State<DenseMatrix, SparseMatrixPolicy, LuDecompositionPolicy>>;

}  // namespace micm

#include "solver_builder.inl"
