// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_info.hpp>
#include <micm/external_model.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cstddef>
#include <functional>
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

    /// @brief Pre-compiled constraint parameter functions (initialized during solver build via SetConstraintFunctions)
    std::vector<std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)>> constraint_param_functions_;

    /// @brief Pre-compiled constraint residual functions (initialized during solver build via SetConstraintFunctions)
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>>
        constraint_forcing_functions_;

    /// @brief Pre-compiled constraint Jacobian functions (initialized during solver build via SetConstraintFunctions)
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>>
        constraint_jacobian_functions_;

    /// @brief External model constraint wrappers
    std::vector<ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> external_constraint_models_;

    /// @brief Runtime count of algebraic variables contributed by external models
    std::size_t external_constraint_count_ = 0;

    /// @brief Pre-compiled external constraint residual functions
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>>
        external_constraint_forcing_functions_;

    /// @brief Pre-compiled external constraint Jacobian functions
    std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>>
        external_constraint_jacobian_functions_;

    /// @brief Pre-compiled external constraint parameter update functions
    std::vector<std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)>>
        external_constraint_param_functions_;

    /// @brief Pre-compiled external constraint parameter initialization functions
    ///        These diagnose constraint parameters from state variables at the start of each Solve()
    std::vector<std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>> external_constraint_init_functions_;

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

    /// @brief Get the number of constraints (built-in + external model)
    std::size_t Size() const
    {
      return constraints_.size() + external_constraint_count_;
    }

    /// @brief Returns species ids whose rows are algebraic when constraints replace state rows
    /// @return Set of variable ids for algebraic rows
    const std::set<std::size_t>& AlgebraicVariableIds() const
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
      std::unordered_map<std::string, int> name_counts;

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

      for (auto& each : constraints_)
      {
        std::visit(
            [&](auto& c)
            {
              for (auto& label : c.parameters_)
                param_names.insert(label);
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
      for (const auto& forcing_fn : constraint_forcing_functions_)
        forcing_fn(state_variables, state_parameters, forcing);
      for (const auto& forcing_fn : external_constraint_forcing_functions_)
        forcing_fn(state_variables, state_parameters, forcing);
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
      for (const auto& jacobian_fn : constraint_jacobian_functions_)
        jacobian_fn(state_variables, state_parameters, jacobian);
      for (const auto& jacobian_fn : external_constraint_jacobian_functions_)
        jacobian_fn(state_variables, state_parameters, jacobian);
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
        const auto& state_variable_indices,   // std::unordered_map<std::string, std::size_t>
        const auto& state_parameter_indices,  // std::unordered_map<std::string, std::size_t>
        SparseMatrixPolicy& jacobian)
    {
      SetConstraintParamIndices(state_parameter_indices);

      constraint_param_functions_.clear();
      constraint_forcing_functions_.clear();
      constraint_jacobian_functions_.clear();

      for (const auto& info : constraint_info_)
      {
        constraint_param_functions_.push_back(
            constraints_[info.index_].template ConstraintParameterFunction<DenseMatrixPolicy>(info));

        constraint_forcing_functions_.push_back(constraints_[info.index_].template ResidualFunction<DenseMatrixPolicy>(
            info, state_variable_indices, state_parameter_indices));

        constraint_jacobian_functions_.push_back(
            constraints_[info.index_].template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
                info,
                state_variable_indices,
                state_parameter_indices,
                jacobian_flat_ids_.begin() + info.jacobian_flat_offset_,
                jacobian));
      }
    }

    /// @brief Returns pre-compiled constraint parameter update functions
    ///        These functions compute temperature-dependent parameters (e.g., K_eq) for each constraint
    ///        Called by solver builder to retrieve functions for UpdateStateParameters pipeline
    /// @return Vector of function objects that take (conditions, state_param) and update constraint parameters
    auto GetUpdateStateParamFunctions()
    {
      return constraint_param_functions_;
    }

    /// @brief Set external model constraint wrappers
    /// @param models Vector of type-erased external model constraint wrappers
    void SetExternalConstraintModels(std::vector<ExternalModelConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>>&& models)
    {
      external_constraint_models_ = std::move(models);
    }

    /// @brief Resolve external model constraints at runtime
    ///
    /// Calls each external model's `algebraic_variable_names_func_()` to determine which
    /// (if any) algebraic variables it contributes. Models returning empty sets are skipped.
    /// Populates `external_constraint_count_` and adds to `algebraic_variable_ids_`.
    ///
    /// @param variable_map Map from species names to state variable indices
    void ResolveExternalConstraints(const std::unordered_map<std::string, std::size_t>& variable_map)
    {
      external_constraint_count_ = 0;
      for (const auto& model : external_constraint_models_)
      {
        auto alg_names = model.algebraic_variable_names_func_();
        for (const auto& name : alg_names)
        {
          auto it = variable_map.find(name);
          if (it == variable_map.end())
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
                "External model constraint targets unknown algebraic species '" + name + "'");
          }
          if (!algebraic_variable_ids_.insert(it->second).second)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_ALGEBRAIC_SPECIES,
                "Multiple constraints map to the same algebraic species row '" + name + "'");
          }
          ++external_constraint_count_;
        }
      }
    }

    /// @brief Returns non-zero Jacobian elements contributed by external model constraints
    /// @param variable_map Map from species names to state variable indices
    /// @return Set of (row, column) index pairs
    std::set<std::pair<std::size_t, std::size_t>> ExternalNonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& variable_map) const
    {
      std::set<std::pair<std::size_t, std::size_t>> ids;
      for (const auto& model : external_constraint_models_)
      {
        auto alg_names = model.algebraic_variable_names_func_();
        if (alg_names.empty())
          continue;
        auto model_ids = model.non_zero_jacobian_elements_func_(variable_map);
        ids.insert(model_ids.begin(), model_ids.end());
        // Ensure diagonal elements exist for algebraic rows
        for (const auto& name : alg_names)
        {
          auto it = variable_map.find(name);
          if (it != variable_map.end())
            ids.insert(std::make_pair(it->second, it->second));
        }
      }
      return ids;
    }

    /// @brief Returns all unique state parameter names from external constraint models
    /// @return Vector of parameter names (duplicates across models are detected and rejected)
    std::vector<std::string> ExternalConstraintParameterNames() const
    {
      std::set<std::string> seen;
      std::vector<std::string> names;
      for (const auto& model : external_constraint_models_)
      {
        auto alg_names = model.algebraic_variable_names_func_();
        if (alg_names.empty())
          continue;
        auto param_names = model.state_parameter_names_func_();
        for (const auto& name : param_names)
        {
          if (!seen.insert(name).second)
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate external constraint parameter name across models: " + name);
          names.push_back(name);
        }
      }
      return names;
    }

    /// @brief Returns pre-compiled external constraint parameter update functions
    auto GetExternalUpdateStateParamFunctions() const
    {
      return external_constraint_param_functions_;
    }

    /// @brief Returns all unique state parameter names that need initialization from state variables
    /// @return Vector of parameter names for state-diagnosed constraint parameters
    std::vector<std::string> ExternalInitializeConstraintParameterNames() const
    {
      std::set<std::string> seen;
      std::vector<std::string> names;
      for (const auto& model : external_constraint_models_)
      {
        auto alg_names = model.algebraic_variable_names_func_();
        if (alg_names.empty())
          continue;
        auto init_names = model.initialize_constraint_parameter_names_func_();
        for (const auto& name : init_names)
        {
          if (!seen.insert(name).second)
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate external initialize constraint parameter name across models: " + name);
          names.push_back(name);
        }
      }
      return names;
    }

    /// @brief Returns pre-compiled external constraint parameter initialization functions
    auto GetExternalInitializeConstraintParamFunctions() const
    {
      return external_constraint_init_functions_;
    }

    /// @brief Initializes constraint parameters from current state variables
    ///        Called at the beginning of each Solve() to diagnose state-dependent constraint constants
    /// @param state_variables Current species concentrations
    /// @param state_parameters State parameters to be updated with diagnosed values
    void InitializeConstraintParameters(const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& state_parameters) const
    {
      for (const auto& init_func : external_constraint_init_functions_)
        init_func(state_variables, state_parameters);
    }

    /// @brief Pre-compiles external constraint residual, Jacobian, and parameter update functions
    ///        Must be called after ResolveExternalConstraints and after Jacobian is built.
    /// @param state_parameter_indices Map from parameter names to state parameter indices
    /// @param state_variable_indices Map from species names to state variable indices
    /// @param jacobian The sparse Jacobian matrix
    void SetExternalModelConstraintFunctions(
        const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const SparseMatrixPolicy& jacobian)
    {
      external_constraint_forcing_functions_.clear();
      external_constraint_jacobian_functions_.clear();
      external_constraint_param_functions_.clear();
      external_constraint_init_functions_.clear();
      for (const auto& model : external_constraint_models_)
      {
        auto alg_names = model.algebraic_variable_names_func_();
        if (alg_names.empty())
          continue;
        external_constraint_param_functions_.push_back(model.update_state_parameters_function_(state_parameter_indices));
        external_constraint_forcing_functions_.push_back(
            model.get_residual_function_(state_parameter_indices, state_variable_indices));
        external_constraint_jacobian_functions_.push_back(
            model.get_jacobian_function_(state_parameter_indices, state_variable_indices, jacobian));
        external_constraint_init_functions_.push_back(
            model.get_initialize_constraint_parameters_function_(state_parameter_indices, state_variable_indices));
      }
    }

   private:
    /// @brief Maps constraint parameter names to their column indices in the state parameter matrix
    ///        Populates constraint_info_[i].state_param_indices_ for each constraint
    ///        Called internally by SetConstraintFunctions after parameter map is finalized
    /// @param state_parameter_indices Map from parameter names to column indices in state_param matrix
    void SetConstraintParamIndices(const auto& state_parameter_indices)  // std::unordered_map<std::string, std::size_t>)
    {
      for (std::size_t i = 0; i < constraints_.size(); ++i)
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