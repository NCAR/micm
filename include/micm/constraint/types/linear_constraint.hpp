// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/micm_exception.hpp>
#include <micm/constraint/constraint_info.hpp>
#include <micm/system/stoich_species.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <system_error>
#include <vector>

namespace micm
{

  /// @brief Constraint for linear relationships: sum(coeff[i] * [species[i]]) = constant
  ///        For example: A + B + C = 1.0 represents a conservation law
  ///        The linear constraint is: G = c1*[A] + c2*[B] + c3*[C] - constant = 0
  class LinearConstraint
  {
   public:
    /// @brief Name of the constraint
    std::string name_;

    /// @brief Names of species this constraint depends on
    std::vector<std::string> species_dependencies_;

    /// @brief Species and their coefficients in the linear sum
    std::vector<StoichSpecies> terms_;

    /// @brief The constant value the linear sum should equal
    double constant_;

   public:
    /// @brief Default constructor
    LinearConstraint() = default;

    /// @brief Construct a linear constraint
    ///        Validates that terms are non-empty
    ///        Builds species_dependencies_ from terms
    /// @param name Constraint identifier
    /// @param terms Vector of StoichSpecies (species, coefficient) in the linear sum
    /// @param constant The value that sum(coeff[i] * [species[i]]) should equal
    LinearConstraint(
        const std::string& name,
        const std::vector<StoichSpecies>& terms,
        double constant)
        : name_(name),
          terms_(terms),
          constant_(constant)
    {
      if (terms_.empty())
      {
        throw MicmException(MicmSeverity::Error, MICM_ERROR_CATEGORY_CONSTRAINT, MICM_CONSTRAINT_ERROR_CODE_EMPTY_REACTANTS, "");
      }
      for (const auto& term : terms_)
      {
        species_dependencies_.push_back(term.species_.name_);
      }
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    ///        For linear constraints, the last species is used in the terms list as the algebraic row target.
    /// @return Species name of the primary algebraic variable
    const std::string& AlgebraicSpecies() const
    {
      return terms_.back().species_.name_;
    }

    /// @brief Create a function to compute the linear constraint residual
    ///        Returns a reusable function object that evaluates:
    ///        G = sum(coeff[i] * [species[i]]) - constant
    /// @param info Constraint information including species indices and row index
    /// @param state_variable_indices Mapping of state variable names to indices
    /// @return Function object that takes (state_variables, forcing) and computes residual
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>
    ResidualFunction(
      const ConstraintInfo& info,
      const auto& state_variable_indices) const
    {
      DenseMatrixPolicy temp_state_variables{1, state_variable_indices.size(), 0.0};

      // Copy data to avoid issues when ConstraintSet is moved
      std::vector<double> coeffs;
      std::vector<std::size_t> species_indices;
      for (std::size_t i = 0; i < this->terms_.size(); ++i)
      {
        coeffs.push_back(this->terms_[i].coefficient_);
        species_indices.push_back(info.state_indices_[i]);
      }
      double constant = this->constant_;
      std::size_t row_idx = info.row_index_;

      return DenseMatrixPolicy::Function(
        [coeffs, species_indices, constant, row_idx](auto&& state, auto&& force)
        {
          // Create a variable for accumulating the linear sum
          auto linear_sum = force.GetRowVariable();
          
          // Initialize linear_sum to 0.0 before accumulation
          state.ForEachRow([](double& sum)
          {
            sum = 0.0;
          }, linear_sum);

          for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const double coeff = coeffs[i];
            const std::size_t species_idx = species_indices[i];
            
            state.ForEachRow([coeff](const double& conc, double& sum)
            {
              sum += coeff * conc;
            }, state.GetConstColumnView(species_idx), linear_sum);
          }
          
          // Forcing term = sum(coeff[i] * [species[i]]) - constant
          state.ForEachRow([constant](const double& sum_val, double& forcing_term)
          {
            forcing_term = sum_val - constant;
          }, linear_sum, force.GetColumnView(row_idx));
        },
        temp_state_variables, 
        temp_state_variables
      );
    }

    /// @brief Create a function to compute Jacobian entries dG/d[species]
    ///        For a linear constraint, the Jacobian is simply the coefficients:
    ///        dG/d[species[i]] = coeff[i]
    /// @param info Constraint information including species indices
    /// @param state_variable_indices Mapping of state variable names to indices
    /// @param jacobian_flat_ids Iterator to this constraint's flat Jacobian indices in shared storage
    /// @param jacobian Sparse matrix to store Jacobian values
    /// @return Function object that takes (state_variables, jacobian) and computes partials
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)>
    JacobianFunction(
      const ConstraintInfo& info,
      const auto& state_variable_indices,
      auto jacobian_flat_ids,
      SparseMatrixPolicy& jacobian) const
    {
      DenseMatrixPolicy temp_state_variables{1, state_variable_indices.size(), 0.0};

      // Pre-compute flat IDs and store them in a vector
      // This avoids iterator issues when the lambda is called multiple times
      std::vector<std::size_t> flat_ids;
      flat_ids.reserve(this->terms_.size());
      for (std::size_t i = 0; i < this->terms_.size(); ++i)
      {
        flat_ids.push_back(*jacobian_flat_ids++);
      }

      // Copy data to avoid issues when ConstraintSet is moved
      std::vector<double> coeffs;
      for (std::size_t i = 0; i < this->terms_.size(); ++i)
      {
        coeffs.push_back(this->terms_[i].coefficient_);
      }

      return SparseMatrixPolicy::Function(
        [coeffs, flat_ids](auto&& state, auto&& jacobian_values)
        {
          // For linear constraints, dG/d[species[i]] = coeff[i]
          // We subtract the coefficient from the Jacobian (matching the SubtractJacobianTerms convention)
          for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const double coeff = coeffs[i];
            
            jacobian_values.ForEachBlock(
              [coeff](double& jac)
              {
                jac -= coeff;
              },
              jacobian_values.GetBlockView(flat_ids[i])
            );
          }
        },
        temp_state_variables,
        jacobian
      );
    }
  };

}  // namespace micm
