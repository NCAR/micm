// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_info.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Manages a collection of algebraic constraints for DAE solvers
  ///        ConstraintSet handles the computation of constraint residuals (forcing terms)
  ///        and Jacobian contributions for a set of constraints. It follows the same
  ///        pattern as ProcessSet for integration with the Rosenbrock solver.
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  class ConstraintSet
  {
   private:
    /// @brief The constraints
    std::vector<Constraint> constraints_;

    /// @brief Information about each constraint for forcing/Jacobian computation
    std::vector<ConstraintInfo> constraint_info_;

    /// @brief Flat list of species indices for each constraint's dependencies
    std::vector<std::size_t> dependency_ids_;

    /// @brief Flat indices into the Jacobian sparse matrix for each constraint's Jacobian entries
    std::vector<std::size_t> jacobian_flat_ids_;

    /// @brief Species variable ids whose ODE rows are replaced by constraints
    std::set<std::size_t> algebraic_variable_ids_;

    /// @brief Pre-compiled constraint residual functions (initialized during solver build via SetConstraintFunctions)
    std::vector<std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>> constraint_forcing_functions_;

    /// @brief Pre-compiled constraint Jacobian functions (initialized during solver build via SetConstraintFunctions)
    std::vector<std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)>> constraint_jacobian_functions_;

   public:
    /// @brief Default constructor
    ConstraintSet() = default;

    /// @brief Construct a ConstraintSet from constraints and variable mapping
    ///        Constraints replace selected species rows in the state/Jacobian (DAE formulation)
    /// @param constraints Vector of constraints
    /// @param variable_map Map from species names to state variable indices
    ConstraintSet(std::vector<Constraint>&& constraints, const std::unordered_map<std::string, std::size_t>& variable_map)
        : constraints_(std::move(constraints))
    {
      // Build constraint info and dependency indices
      std::size_t dependency_offset = 0;
      for (std::size_t i = 0; i < constraints_.size(); ++i)
      {
        const auto& constraint = constraints_[i];

        ConstraintInfo info;
        info.index_ = i;

        const auto& algebraic_species = constraint.AlgebraicSpecies();
        auto row_it = variable_map.find(algebraic_species);
        if (row_it == variable_map.end())
        {
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
              "Constraint '" + constraint.GetName() + "' targets unknown algebraic species '" + algebraic_species + "'");
        }
        info.row_index_ = row_it->second;

        if (!algebraic_variable_ids_.insert(info.row_index_).second)
        {
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_INVALID_STOICHIOMETRY,
              "Multiple constraints map to the same algebraic species row '" + algebraic_species + "'");
        }

        info.number_of_dependencies_ = constraint.NumberOfDependencies();
        info.dependency_offset_ = dependency_offset;
        info.jacobian_flat_offset_ = 0;  // Set later in SetJacobianFlatIds

        // Map species dependencies to variable indices
        for (const auto& species_name : constraint.SpeciesDependencies())
        {
          auto it = variable_map.find(species_name);
          if (it == variable_map.end())
          {
            throw MicmException(
                MicmSeverity::Error,
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
                "Constraint '" + constraint.GetName() + "' depends on unknown species '" + species_name + "'");
          }
          dependency_ids_.push_back(it->second);
          info.state_indices_.push_back(it->second);  // Also store in ConstraintInfo
        }

        dependency_offset += info.number_of_dependencies_;
        constraint_info_.push_back(info);
      }
    }

    /// @brief Move constructor - default implementation
    ConstraintSet(ConstraintSet&& other) noexcept = default;

    /// @brief Move assignment operator
    ConstraintSet& operator=(ConstraintSet&& other) noexcept = default;

    /// @brief Copy constructor
    ConstraintSet(const ConstraintSet&) = default;

    /// @brief Copy assignment
    ConstraintSet& operator=(const ConstraintSet&) = default;

    /// @brief Get the number of constraints
    std::size_t Size() const
    {
      return constraints_.size();
    }

    /// @brief Returns species ids whose rows are algebraic when constraints replace state rows
    /// @return Set of variable ids for algebraic rows
    const std::set<std::size_t>& AlgebraicVariableIds() const
    {
      return algebraic_variable_ids_;
    }

    /// @brief Add constraint residuals to forcing vector (constraint rows)
    ///        For each constraint G_i, writes or adds G_i(x) to forcing[constraint_row]
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param forcing Forcing terms (grid cell, state variable) - constraint rows will be modified
    void AddForcingTerms(const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing) const
    {
      for (const auto& forcing_fn : constraint_forcing_functions_)
        forcing_fn(state_variables, forcing);
    }

    /// @brief Subtract constraint Jacobian terms from Jacobian matrix
    ///        For each constraint G_i, subtracts dG_i/dx_j from jacobian[constraint_row, j]
    ///        (Subtraction matches the convention used by ProcessSet)
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param jacobian Sparse Jacobian matrix (grid cell, row, column)
    void SubtractJacobianTerms(const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) const
    {
      for (const auto& jacobian_fn : constraint_jacobian_functions_)
        jacobian_fn(state_variables, jacobian);
    }

    /// @brief Returns positions of all non-zero Jacobian elements for constraint rows
    /// @return Set of (row, column) index pairs
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const
    {
      std::set<std::pair<std::size_t, std::size_t>> ids;

      auto dep_id = dependency_ids_.begin();
      for (const auto& info : constraint_info_)
      {
        // Ensure the diagonal element exists for the constraint row (required by AlphaMinusJacobian and LU decomposition)
        ids.insert(std::make_pair(info.row_index_, info.row_index_));
        // Each constraint contributes Jacobian entries at (constraint_row, dependency_column)
        for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
        {
          ids.insert(std::make_pair(info.row_index_, dep_id[i]));
        }
        dep_id += info.number_of_dependencies_;
      }

      return ids;
    }

    /// @brief Computes and stores flat indices for Jacobian elements
    /// @param matrix The sparse Jacobian matrix
    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix)
    {
      jacobian_flat_ids_.clear();

      std::size_t flat_offset = 0;
      for (auto& info : constraint_info_)
      {
        info.jacobian_flat_offset_ = flat_offset;

        // Store flat indices for each dependency of this constraint
        const std::size_t* dep_id = dependency_ids_.data() + info.dependency_offset_;
        for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
        {
          jacobian_flat_ids_.push_back(matrix.VectorIndex(0, info.row_index_, dep_id[i]));
        }

        flat_offset += info.number_of_dependencies_;
      }
    }

    /// @brief Pre-compiles constraint residual and Jacobian functions for efficient evaluation
    ///        Creates reusable function objects from each constraint's ResidualFunction and JacobianFunction.
    ///        Must be called after SetJacobianFlatIds and before solver execution.
    /// @param state_variable_indices Map from species names to state variable indices
    /// @param jacobian The sparse Jacobian matrix (used for function template instantiation)
    void SetConstraintFunctions(
        const auto& state_variable_indices,  // acts like std::unordered_map<std::string, std::size_t>
        SparseMatrixPolicy& jacobian)
    {
      constraint_forcing_functions_.clear();
      constraint_jacobian_functions_.clear();
      for (const auto& info : constraint_info_)
      {
        constraint_forcing_functions_.push_back(
            constraints_[info.index_].template ResidualFunction<DenseMatrixPolicy>(info, state_variable_indices));

        constraint_jacobian_functions_.push_back(
            constraints_[info.index_].template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
                info, state_variable_indices, jacobian_flat_ids_.begin() + info.jacobian_flat_offset_, jacobian));
      }
    }
  };

}  // namespace micm