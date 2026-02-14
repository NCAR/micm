// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/solver.hpp>
#include <micm/solver/solver_error.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <system_error>
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
    bool ignore_unused_species_ = true;
    bool reorder_state_ = true;
    bool valid_system_ = false;

   public:
    SolverBuilder() = delete;

    SolverBuilder(const SolverParametersPolicy& options)
        : options_(options)
    {
    }

    virtual ~SolverBuilder() = default;

    /// @brief Set the chemical system
    /// @param system The chemical system
    /// @return Updated SolverBuilder
    SolverBuilder& SetSystem(const System& system);

    /// @brief Set the reactions
    /// @param reactions The reactions
    /// @return Updated SolverBuilder
    SolverBuilder& SetReactions(const std::vector<Process>& reactions);

    /// @brief Add processes from an external model
    /// @param model The external model
    /// @return Updated SolverBuilder
    template<class ExternalModel>
    SolverBuilder& AddExternalModelProcesses(ExternalModel&& model);

    /// @brief Set whether to ignore unused species
    /// @param ignore_unused_species True if unused species should be ignored
    /// @return Updated SolverBuilder
    SolverBuilder& SetIgnoreUnusedSpecies(bool ignore_unused_species);

    /// @brief Set whether to reorder the state to optimize the LU decomposition
    /// @param reorder_state True if the state should be reordered
    /// @return Updated SolverBuilder
    SolverBuilder& SetReorderState(bool reorder_state);

    /// @brief Creates an instance of Solver with a properly configured ODE solver
    /// @return An instance of Solver
    auto Build() const;

   protected:
    /// @brief Checks for unused species
    /// @param rates The rates policy instance containing information about processes
    /// @throws std::system_error if an unused species is found
    void UnusedSpeciesCheck(const RatesPolicy& rates) const;

    /// @brief Gets a map of species to their index
    /// @return The species map
    std::unordered_map<std::string, std::size_t> GetSpeciesMap() const;

    /// @brief Sets the absolute tolerances per species
    /// @param parameters
    /// @param species_map
    void SetAbsoluteTolerances(
        std::vector<double>& tolerances,
        const std::unordered_map<std::string, std::size_t>& species_map) const;

    /// @brief Returns the labels of the custom parameters
    /// @return The labels of the custom parameters
    std::unordered_map<std::string, std::size_t> GetCustomParameterMap() const;
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