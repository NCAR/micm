// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_info.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/types.hpp>

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
    template<class U>
    using Vector = typename DenseMatrixPolicy::template VectorType<U>;
    template<class U>
    using VectorView = typename DenseMatrixPolicy::template VectorType<U>::ConstViewType;

    /// @brief The constraints
    std::vector<Constraint<DenseMatrixPolicy, SparseMatrixPolicy>> constraints_;

    /// @brief Information about each constraint for forcing/Jacobian computation
    std::vector<ConstraintInfo> constraint_info_;

    /// @brief Flat list of species indices for each constraint's dependencies
    std::vector<Index> dependency_ids_;

    /// @brief Flat indices into the Jacobian sparse matrix for each constraint's Jacobian entries
    std::vector<Index> jacobian_flat_ids_;

    /// @brief Species variable ids whose ODE rows are replaced by constraints
    std::set<Index> algebraic_variable_ids_;

    /// @brief Pre-compiled function to set algebraic variable error estimates
    ///        For algebraic variables, the embedded error formula produces near-zero Yerror (because M_ii = 0).
    ///        Instead, we use the step change (Ynew[a] - Y[a]) as the error estimate for algebraic variables.
    ///        This captures the full change imposed by constraint enforcement and allows the solver to
    ///        reject steps where algebraic variables change too much relative to their tolerances.
    Vector<Index> alg_ids_data_;
    VectorView<Index> alg_ids_view_;

   public:
    /// @brief Default constructor
    ConstraintSet() = default;

    /// @brief Construct a ConstraintSet from constraints and variable mapping
    ///        Constraints replace selected species rows in the state/Jacobian (DAE formulation)
    /// @param constraints Vector of constraints
    /// @param variable_map Map from species names to state variable indices
    ConstraintSet(
        std::vector<Constraint<DenseMatrixPolicy, SparseMatrixPolicy>>&& constraints,
        const std::unordered_map<std::string, Index>& variable_map)
        : constraints_(std::move(constraints))
    {
      // Build constraint info and dependency indices
      Index dependency_offset = 0;
      for (Index i = 0; i < constraints_.size(); ++i)
      {
        const auto& constraint = constraints_[i];

        ConstraintInfo info;
        info.index_ = i;

        const auto& algebraic_species = constraint.AlgebraicSpecies();
        auto row_it = variable_map.find(algebraic_species);
        if (row_it == variable_map.end())
        {
          throw MicmException(
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
              "Constraint '" + constraint.GetName() + "' targets unknown algebraic species '" + algebraic_species + "'");
        }
        info.row_index_ = row_it->second;

        if (!algebraic_variable_ids_.insert(info.row_index_).second)
        {
          throw MicmException(
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
    Index Size() const
    {
      return constraints_.size();
    }

    /// @brief Returns species ids whose rows are algebraic when constraints replace state rows
    /// @return Set of variable ids for algebraic rows
    const std::set<Index>& AlgebraicVariableIds() const
    {
      return algebraic_variable_ids_;
    }

    /// @brief Deduplicates parameter names across all constraints in the set
    ///        Ensures all constraint parameters have globally unique names by appending
    ///        numeric suffixes (_1, _2, etc.) to duplicates.
    ///        This should be called immediately after construction so that parameter
    ///        names are finalized before the solver builder creates the parameter map.
    ///        This logic is not part of the constructor because it mutates the constraint
    ///        parameters, which is considered beyond the scope of construction.
    void SetUniqueParameterNames()
    {
      std::unordered_set<std::string> used_names;
      std::unordered_map<std::string, Index> name_counts;

      for (auto& each : constraints_)
      {
        std::visit(
            [&](auto& c)
            {
              for (auto& label : c.parameters_)
              {
                const std::string original = label;

                if (used_names.count(label) > 0)
                {
                  auto& count = name_counts[original];
                  do
                  {
                    count++;
                    label = original + "_" + std::to_string(count);
                  } while (used_names.count(label) > 0);
                }
                used_names.insert(label);
              }
            },
            each.constraint_);
      }
    }

    /// @brief Returns all unique parameter names from all constraints in the set
    /// @return Set of parameter names
    std::unordered_set<std::string> GetParameterNames() const
    {
      std::unordered_set<std::string> param_names;

      for (const auto& each : constraints_)
      {
        std::visit(
            [&](auto& c)
            {
              for (auto& label : c.parameters_)
              {
                param_names.insert(label);
              }
            },
            each.constraint_);
      }

      return param_names;
    }

    /// @brief Add constraint residuals to forcing vector (constraint rows)
    ///        For each constraint G_i, writes or adds G_i(x) to forcing[constraint_row]
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param state_parameters Current state parameters (grid cell, parameter) - e.g., temperature-dependent K_eq values
    /// @param forcing Forcing terms (grid cell, state variable) - constraint rows will be modified
    void AddForcingTerms(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& forcing) const
    {
      for (const auto& info : constraint_info_)
      {
        constraints_[info.index_].AddResidual(info, state_variables, state_parameters, forcing);
      }
    }

    /// @brief Subtract constraint Jacobian terms from Jacobian matrix
    ///        For each constraint G_i, subtracts dG_i/dx_j from jacobian[constraint_row, j]
    ///        (Subtraction matches the convention used by ProcessSet)
    /// @param state_variables Current species concentrations (grid cell, species)
    /// @param state_parameters Current state parameters (grid cell, parameter) - e.g., temperature-dependent K_eq values
    /// @param jacobian Sparse Jacobian matrix (grid cell, row, column)
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      for (const auto& info : constraint_info_)
      {
        constraints_[info.index_].SubtractJacobian(info, state_variables, state_parameters, jacobian);
      }
    }

    /// @brief Set algebraic variable error estimates using step changes
    ///        For each algebraic variable a: Yerror[a] = Ynew[a] - Y[a]
    /// @param Yerror Error vector — algebraic entries are overwritten with step changes
    /// @param Y State at beginning of step
    /// @param Ynew Proposed state at end of step (after constraint enforcement)
    void SetAlgebraicErrors(DenseMatrixPolicy& Yerror, const DenseMatrixPolicy& Y, const DenseMatrixPolicy& Ynew) const
    {
      if (algebraic_variable_ids_.empty())
      {
        return;
      }

      const auto& alg_ids = alg_ids_view_;

      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& yerr,
              const typename DenseMatrixPolicy::ConstViewType& y,
              const typename DenseMatrixPolicy::ConstViewType& ynew) {
            for (const auto& col : alg_ids)
            {
              yerr.ForEachRow(
                  [](const Real& ynew_a, const Real& y_a, Real& err_a) { err_a = ynew_a - y_a; },
                  ynew.GetConstColumnView(col),
                  y.GetConstColumnView(col),
                  yerr.GetColumnView(col));
            }
          },
          Yerror,
          Y,
          Ynew)(Yerror, Y, Ynew);
    }

    /// @brief Returns positions of all non-zero Jacobian elements for constraint rows
    /// @return Set of (row, column) index pairs
    std::set<std::pair<Index, Index>> NonZeroJacobianElements() const
    {
      std::set<std::pair<Index, Index>> ids;

      auto dep_id = dependency_ids_.begin();
      for (const auto& info : constraint_info_)
      {
        // Ensure the diagonal element exists for the constraint row (required by AlphaMinusJacobian and LU decomposition)
        ids.insert(std::make_pair(info.row_index_, info.row_index_));
        // Each constraint contributes Jacobian entries at (constraint_row, dependency_column)
        for (Index i = 0; i < info.number_of_dependencies_; ++i)
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
    void SetJacobianFlatIds(const SparseMatrix<Real, OrderingPolicy>& matrix)
    {
      jacobian_flat_ids_.clear();

      Index flat_offset = 0;
      for (auto& info : constraint_info_)
      {
        info.jacobian_flat_offset_ = flat_offset;

        // Store flat indices for each dependency of this constraint
        const Index* dep_id = dependency_ids_.data() + info.dependency_offset_;
        for (Index i = 0; i < info.number_of_dependencies_; ++i)
        {
          jacobian_flat_ids_.push_back(matrix.VectorIndex(0, info.row_index_, dep_id[i]));
        }

        flat_offset += info.number_of_dependencies_;
      }
    }

    /// @brief Sets up constraint indices and Jacobian metadata for direct-call execution.
    ///        Must be called after SetJacobianFlatIds and before solver execution.
    /// @param state_parameter_indices Map from parameter names to state parameter indices
    void SetConstraintFunctions(const auto& state_parameter_indices)
    {
      SetConstraintParamIndices(state_parameter_indices);

      for (const auto& info : constraint_info_)
      {
        auto flat_ids = jacobian_flat_ids_.begin() + info.jacobian_flat_offset_;
        constraints_[info.index_].SetStateIndices(info, flat_ids);
      }

      BuildAlgebraicErrorFunction();
    }

    /// @brief Apply constraint parameter updates for all grid cells (e.g., temperature-dependent K_eq).
    ///        Called directly from the solver's UpdateStateParameters pipeline.
    /// @param conditions Per-grid-cell atmospheric conditions
    /// @param state_param State parameter matrix to update
    void UpdateStateParameters(
        const typename DenseMatrixPolicy::template VectorType<Conditions>& conditions,
        DenseMatrixPolicy& state_param) const
    {
      for (const auto& info : constraint_info_)
      {
        constraints_[info.index_].ApplyConstraintParameter(info, conditions, state_param);
      }
    }

    /// @brief Extend the algebraic variable set with rows contributed by external models.
    ///        Called by SolverBuilder after collecting external algebraic-variable IDs so that
    ///        Yerror handling covers external algebraic rows as well.
    void AddExternalAlgebraicVariableIds(const std::set<Index>& ids)
    {
      algebraic_variable_ids_.insert(ids.begin(), ids.end());
    }

    /// @brief Rebuilds the algebraic-variable id view used by SetAlgebraicErrors.
    ///        Call after any external algebraic ids have been merged in.
    void FinalizeAlgebraicErrorFunction()
    {
      if (algebraic_variable_ids_.empty())
      {
        return;
      }
      std::vector<Index> alg_ids_temp(algebraic_variable_ids_.begin(), algebraic_variable_ids_.end());
      alg_ids_data_ = alg_ids_temp;
      alg_ids_data_.CopyToDevice();
      alg_ids_view_ = std::as_const(alg_ids_data_).GetView();
    }

   private:
    void BuildAlgebraicErrorFunction()
    {
      FinalizeAlgebraicErrorFunction();
    }

    /// @brief Maps constraint parameter names to their column indices in the state parameter matrix
    ///        Populates constraint_info_[i].state_param_indices_ for each constraint
    ///        Called internally by SetConstraintFunctions after parameter map is finalized
    /// @param state_parameter_indices Map from parameter names to column indices in state_param matrix
    void SetConstraintParamIndices(const auto& state_parameter_indices)  // std::unordered_map<std::string, Index>)
    {
      for (Index i = 0; i < constraints_.size(); ++i)
      {
        for (const auto& name : constraints_[i].GetParameterNames())
        {
          auto it = state_parameter_indices.find(name);
          if (it == state_parameter_indices.end())
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
                "Constraint '" + constraints_[i].GetName() + "' depends on unknown parameter '" + name + "'");
          }
          constraint_info_[i].state_param_indices_.push_back(it->second);  // store in ConstraintInfo
        }
      }
    }
  };

}  // namespace micm