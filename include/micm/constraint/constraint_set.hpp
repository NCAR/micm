// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_error.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cstddef>
#include <unordered_map>
#include <memory>
#include <set>
#include <string>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Manages a collection of algebraic constraints for DAE solvers
  ///        ConstraintSet handles the computation of constraint residuals (forcing terms)
  ///        and Jacobian contributions for a set of constraints. It follows the same
  ///        pattern as ProcessSet for integration with the Rosenbrock solver.
  class ConstraintSet
  {
   protected:
    /// @brief Information for each constraint
    struct ConstraintInfo
    {
      std::size_t constraint_index_;       // Index in constraints_ vector
      std::size_t constraint_row_;         // Row in the forcing/Jacobian
      std::size_t number_of_dependencies_; // Number of species this constraint depends on
      std::size_t dependency_offset_;      // Starting offset in dependency_ids_
      std::size_t jacobian_flat_offset_;   // Starting offset in jacobian_flat_ids_
    };

    /// @brief The constraints
    std::vector<Constraint> constraints_;

    /// @brief Flat list of species indices for each constraint's dependencies
    std::vector<std::size_t> dependency_ids_;

    /// @brief Information about each constraint for forcing/Jacobian computation
    std::vector<ConstraintInfo> constraint_info_;

    /// @brief Flat indices into the Jacobian sparse matrix for each constraint's Jacobian entries
    std::vector<std::size_t> jacobian_flat_ids_;

    /// @brief Row offset for constraints in the extended state (= number of species)
    std::size_t constraint_row_offset_{ 0 };

    /// @brief If true, constraints replace selected species rows instead of augmenting extra rows
    bool replace_state_rows_{ false };

    /// @brief Species variable ids whose ODE rows are replaced by constraints
    std::set<std::size_t> algebraic_variable_ids_;

    /// @brief Maximum number of dependencies across all constraints (for buffer allocation)
    std::size_t max_dependencies_{ 0 };

   public:
    /// @brief Default constructor
    ConstraintSet() = default;

    /// @brief Construct a ConstraintSet from constraints and variable mapping
    /// @param constraints Vector of constraints
    /// @param variable_map Map from species names to state variable indices
    /// @param constraint_row_offset Row offset for augmented-rows mode (= number of species)
    /// @param replace_state_rows If true, constraints replace selected species rows in the state/Jacobian
    ConstraintSet(
        std::vector<Constraint>&& constraints,
        const std::unordered_map<std::string, std::size_t>& variable_map,
        std::size_t constraint_row_offset,
        bool replace_state_rows = false);

    /// @brief Move constructor
    ConstraintSet(ConstraintSet&& other) noexcept = default;

    /// @brief Move assignment
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

    /// @brief Returns positions of all non-zero Jacobian elements for constraint rows
    /// @return Set of (row, column) index pairs
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;

    /// @brief Returns species ids whose rows are algebraic when constraints replace state rows
    /// @return Set of variable ids for algebraic rows
    const std::set<std::size_t>& AlgebraicVariableIds() const
    {
      return algebraic_variable_ids_;
    }

    /// @brief Computes and stores flat indices for Jacobian elements
    /// @param matrix The sparse Jacobian matrix
    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);

    /// @brief Add constraint residuals to forcing vector (constraint rows)
    ///        For each constraint G_i, writes or adds G_i(x) to forcing[constraint_row]
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param forcing Forcing terms (grid cell, state variable) - constraint rows will be modified
    template<typename DenseMatrixPolicy>
    void AddForcingTerms(const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing) const;

    /// @brief Subtract constraint Jacobian terms from Jacobian matrix
    ///        For each constraint G_i, subtracts dG_i/dx_j from jacobian[constraint_row, j]
    ///        (Subtraction matches the convention used by ProcessSet)
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param jacobian Sparse Jacobian matrix (grid cell, row, column)
    template<class DenseMatrixPolicy, class SparseMatrixPolicy>
    void SubtractJacobianTerms(const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian) const;
  };

  inline ConstraintSet::ConstraintSet(
      std::vector<Constraint>&& constraints,
      const std::unordered_map<std::string, std::size_t>& variable_map,
      std::size_t constraint_row_offset,
      bool replace_state_rows)
      : constraints_(std::move(constraints)),
        constraint_row_offset_(constraint_row_offset),
        replace_state_rows_(replace_state_rows)
  {
    // Build constraint info and dependency indices
    std::size_t dependency_offset = 0;
    for (std::size_t i = 0; i < constraints_.size(); ++i)
    {
      const auto& constraint = constraints_[i];

      ConstraintInfo info;
      info.constraint_index_ = i;
      if (replace_state_rows_)
      {
        const auto algebraic_species = constraint.GetAlgebraicSpecies();
        auto row_it = variable_map.find(algebraic_species);
        if (row_it == variable_map.end())
        {
          throw std::system_error(
              make_error_code(MicmConstraintErrc::UnknownSpecies),
              "Constraint '" + constraint.GetName() + "' targets unknown algebraic species '" + algebraic_species + "'");
        }
        info.constraint_row_ = row_it->second;

        if (!algebraic_variable_ids_.insert(info.constraint_row_).second)
        {
          throw std::runtime_error(
              "Multiple constraints map to the same algebraic species row '" + algebraic_species + "'");
        }
      }
      else
      {
        info.constraint_row_ = constraint_row_offset_ + i;
      }
      info.number_of_dependencies_ = constraint.NumberOfDependencies();
      info.dependency_offset_ = dependency_offset;
      info.jacobian_flat_offset_ = 0;  // Set later in SetJacobianFlatIds

      // Track maximum dependencies for buffer allocation
      if (info.number_of_dependencies_ > max_dependencies_)
      {
        max_dependencies_ = info.number_of_dependencies_;
      }

      // Map species dependencies to variable indices
      for (const auto& species_name : constraint.GetSpeciesDependencies())
      {
        auto it = variable_map.find(species_name);
        if (it == variable_map.end())
        {
          throw std::system_error(
              make_error_code(MicmConstraintErrc::UnknownSpecies),
              "Constraint '" + constraint.GetName() + "' depends on unknown species '" + species_name + "'");
        }
        dependency_ids_.push_back(it->second);
      }

      dependency_offset += info.number_of_dependencies_;
      constraint_info_.push_back(info);
    }
  }

  inline std::set<std::pair<std::size_t, std::size_t>> ConstraintSet::NonZeroJacobianElements() const
  {
    std::set<std::pair<std::size_t, std::size_t>> ids;

    auto dep_id = dependency_ids_.begin();
    for (const auto& info : constraint_info_)
    {
      // Ensure the diagonal element exists for the constraint row (required by AlphaMinusJacobian and LU decomposition)
      if (replace_state_rows_)
      {
        ids.insert(std::make_pair(info.constraint_row_, info.constraint_row_));
      }
      // Each constraint contributes Jacobian entries at (constraint_row, dependency_column)
      for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
      {
        ids.insert(std::make_pair(info.constraint_row_, dep_id[i]));
      }
      dep_id += info.number_of_dependencies_;
    }

    return ids;
  }

  template<typename OrderingPolicy>
  inline void ConstraintSet::SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix)
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
        jacobian_flat_ids_.push_back(matrix.VectorIndex(0, info.constraint_row_, dep_id[i]));
      }

      flat_offset += info.number_of_dependencies_;
    }
  }

  template<typename DenseMatrixPolicy>
  inline void ConstraintSet::AddForcingTerms(
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    if (constraints_.empty())
      return;

    if constexpr (VectorizableDense<DenseMatrixPolicy>)
    {
      // Vectorized dense layouts are not row-contiguous, so build a contiguous row buffer.
      std::vector<double> concentrations(state_variables.NumColumns());
      for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
      {
        for (std::size_t i_var = 0; i_var < state_variables.NumColumns(); ++i_var)
        {
          concentrations[i_var] = state_variables[i_cell][i_var];
        }

        auto cell_forcing = forcing[i_cell];
        for (const auto& info : constraint_info_)
        {
          const std::size_t* indices = dependency_ids_.data() + info.dependency_offset_;
          double residual = constraints_[info.constraint_index_].Residual(concentrations.data(), indices);
          if (replace_state_rows_)
          {
            cell_forcing[info.constraint_row_] = residual;
          }
          else
          {
            cell_forcing[info.constraint_row_] += residual;
          }
        }
      }
    }
    else
    {
      // Loop over grid cells
      for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
      {
        auto cell_state = state_variables[i_cell];
        auto cell_forcing = forcing[i_cell];

        // Get pointer to concentration data for this cell
        const double* concentrations = &cell_state[0];

        for (const auto& info : constraint_info_)
        {
          // Get pointer to indices for this constraint
          const std::size_t* indices = dependency_ids_.data() + info.dependency_offset_;

          // Evaluate constraint residual and add to forcing
          double residual = constraints_[info.constraint_index_].Residual(concentrations, indices);
          if (replace_state_rows_)
          {
            cell_forcing[info.constraint_row_] = residual;
          }
          else
          {
            cell_forcing[info.constraint_row_] += residual;
          }
        }
      }
    }
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline void ConstraintSet::SubtractJacobianTerms(
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    if (constraints_.empty())
      return;

    // Allocate reusable buffer for constraint Jacobian values
    std::vector<double> jac_buffer(max_dependencies_);

    [[maybe_unused]] std::vector<double> concentrations;
    if constexpr (VectorizableDense<DenseMatrixPolicy>)
    {
      concentrations.resize(state_variables.NumColumns());
    }

    [[maybe_unused]] auto cell_jacobian = jacobian.AsVector().begin();

    // Loop over grid cells
    for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
    {
      const double* concentration_ptr;
      if constexpr (VectorizableDense<DenseMatrixPolicy>)
      {
        // Vectorized dense layouts are not row-contiguous, so build a contiguous row buffer.
        for (std::size_t i_var = 0; i_var < state_variables.NumColumns(); ++i_var)
        {
          concentrations[i_var] = state_variables[i_cell][i_var];
        }
        concentration_ptr = concentrations.data();
      }
      else
      {
        auto cell_state = state_variables[i_cell];
        concentration_ptr = &cell_state[0];
      }

      for (const auto& info : constraint_info_)
      {
        const std::size_t* indices = dependency_ids_.data() + info.dependency_offset_;
        constraints_[info.constraint_index_].Jacobian(concentration_ptr, indices, jac_buffer.data());

        // In replace mode, the Jacobian row is already zeroed (Fill(0) + ProcessSet skips algebraic rows),
        // so -= is safe and correctly accumulates when a species appears in multiple dependency slots.
        if constexpr (VectorizableSparse<SparseMatrixPolicy>)
        {
          const std::size_t* dep_id = dependency_ids_.data() + info.dependency_offset_;
          for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
          {
            jacobian[i_cell][info.constraint_row_][dep_id[i]] -= jac_buffer[i];
          }
        }
        else
        {
          const std::size_t* flat_ids = jacobian_flat_ids_.data() + info.jacobian_flat_offset_;
          for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
          {
            cell_jacobian[flat_ids[i]] -= jac_buffer[i];
          }
        }
      }

      if constexpr (!VectorizableSparse<SparseMatrixPolicy>)
      {
        // Advance to next grid cell's Jacobian block
        cell_jacobian += jacobian.FlatBlockSize();
      }
    }
  }

}  // namespace micm
