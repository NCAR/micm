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
#include <micm/util/types.hpp>
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
    std::vector<Constraint> constraints_;

    std::vector<ExternalModelSystem> external_systems_;
    std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_process_sets_;
    std::vector<ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_constraints_;

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

    SolverBuilder(const SolverBuilder&) = default;
    SolverBuilder& operator=(const SolverBuilder&) = default;
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

    /// @brief Add an external model (state variables, processes, and/or constraints)
    ///
    /// If the model satisfies HasState, its state variables and parameters are registered with the solver.
    /// The model must satisfy at least one of HasProcesses (process wrappers are created) or
    /// HasConstraints (constraint wrappers are created).
    /// @param model The external model (taken by value; caller decides whether to copy or move)
    /// @return Updated SolverBuilder
    template<class ExternalModel>
    SolverBuilder& AddExternalModel(ExternalModel model)
    {
      static_assert(
          HasProcesses<ExternalModel> || HasConstraints<ExternalModel>,
          "External model passed to AddExternalModel() must satisfy at least HasProcesses or HasConstraints");

      if constexpr (HasState<ExternalModel>)
      {
        external_systems_.emplace_back(ExternalModelSystem{ model });
      }

      if constexpr (HasProcesses<ExternalModel> && HasConstraints<ExternalModel>)
      {
        external_process_sets_.emplace_back(ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ model });
        external_constraints_.emplace_back(
            ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::move(model) });
      }
      else if constexpr (HasProcesses<ExternalModel>)
      {
        external_process_sets_.emplace_back(
            ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::move(model) });
      }
      else
      {
        external_constraints_.emplace_back(
            ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::move(model) });
      }
      return *this;
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
    void SetAbsoluteTolerances(
        std::vector<Real>& tolerances,
        const std::unordered_map<std::string, Index>& species_map) const;
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
      class DenseMatrixPolicy = Matrix<Real>,
      class SparseMatrixPolicy = SparseMatrix<Real, SparseMatrixStandardOrdering>,
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
      class DenseMatrix = Matrix<Real>,
      class SparseMatrixPolicy = SparseMatrix<Real, SparseMatrixStandardOrdering>,
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
