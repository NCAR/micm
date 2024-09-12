// Copyright (C) 2023-2024 National Center for Atmospheric Research
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
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <system_error>

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
      class LinearSolverPolicy,
      class StatePolicy>
  class SolverBuilder
  {
   public:
    using DenseMatrixPolicyType = DenseMatrixPolicy;
    using SparseMatrixPolicyType = SparseMatrixPolicy;

   protected:
    SolverParametersPolicy options_;
    System system_;
    std::size_t number_of_grid_cells_ = 1;
    std::vector<Process> reactions_;
    bool ignore_unused_species_ = true;
    bool reorder_state_ = true;
    bool valid_system_ = false;
    bool valid_reactions_ = false;

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

    /// @brief Set the number of grid cells
    /// @param number_of_grid_cells The number of grid cells
    /// @return Updated SolverBuilder
    SolverBuilder& SetNumberOfGridCells(int number_of_grid_cells);

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
    /// @throws std::system_error if an unused species is found
    void UnusedSpeciesCheck() const;

    /// @brief Gets a map of species to their index
    /// @return The species map
    std::map<std::string, std::size_t> GetSpeciesMap() const;

    /// @brief Sets the absolute tolerances per species
    /// @param parameters
    /// @param species_map
    void SetAbsoluteTolerances(std::vector<double>& tolerances, const std::map<std::string, std::size_t>& species_map) const;

    /// @brief Returns the labels of the custom parameters
    /// @return The labels of the custom parameters
    std::vector<std::string> GetCustomParameterLabels() const;
  };

  
  /// @brief Builder of CPU-based general solvers
  /// @tparam SolverParametersPolicy Parameters for the ODE solver
  /// @tparam DenseMatrixPolicy Policy for dense matrices
  /// @tparam SparseMatrixPolicy Policy for sparse matrices
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy = Matrix<double>,
      class SparseMatrixPolicy = SparseMatrix<double, SparseMatrixStandardOrdering>>
  using CpuSolverBuilder = SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      ProcessSet,
      LinearSolver<SparseMatrixPolicy, LuDecomposition>,
      State<typename SolverParametersPolicy::template TemporaryVariablesType<DenseMatrixPolicy>, DenseMatrixPolicy, SparseMatrixPolicy>>;
}  // namespace micm

#include "solver_builder.inl"