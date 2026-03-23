// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/micm_exception.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <system_error>
#include <vector>

namespace micm
{
  /// @brief Constraint for chemical equilibrium: K_eq = [products]^stoich / [reactants]^stoich
  ///        For a reversible reaction: aA + bB <-> cC + dD
  ///        The equilibrium constraint is: G = K_eq * [A]^a * [B]^b - [C]^c * [D]^d = 0
  ///        This can also be written in terms of forward/backward rate constants:
  ///        G = k_f * [A]^a * [B]^b - k_b * [C]^c * [D]^d = 0
  ///        where K_eq = k_f / k_b
  class EquilibriumConstraint
  {
   public:
    /// @brief Name of the constraint (for identification)
    std::string name_;

    /// @brief Names of species this constraint depends on
    std::vector<std::string> species_dependencies_;

    /// @brief Reactant species and their stoichiometric coefficients
    std::vector<StoichSpecies> reactants_;

    /// @brief Product species and their stoichiometric coefficients
    std::vector<StoichSpecies> products_;

    /// @brief Equilibrium constant K_eq = k_forward / k_backward
    double equilibrium_constant_;

   private:
    /// @brief Indices into the reactants_ vector for each species dependency
    std::vector<std::size_t> reactant_dependency_indices_;

    /// @brief Indices into the products_ vector for each species dependency
    std::vector<std::size_t> product_dependency_indices_;

   public:
    /// @brief Default constructor
    EquilibriumConstraint() = default;

    /// @brief Construct an equilibrium constraint
    ///        Validates that equilibrium constraint > 0
    ///        Builds species_dependencies_ by concatenating reactants then products
    ///        Stores index mappings for efficient Jacobian computation
    /// @param name Constraint identifier
    /// @param reactants Vector of StoichSpecies (species, stoichiometry) for reactants
    /// @param products Vector of StoichSpecies (species, stoichiometry) for products
    /// @param equilibrium_constant K_eq = [products]/[reactants] at equilibrium
    EquilibriumConstraint(
        const std::string& name,
        const std::vector<StoichSpecies>& reactants,
        const std::vector<StoichSpecies>& products,
        double equilibrium_constant)
        : name_(name),
          reactants_(reactants),
          products_(products),
          equilibrium_constant_(equilibrium_constant)
    {
      if (reactants_.empty())
      {
        throw MicmException(
            MicmSeverity::Error,
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_EMPTY_REACTANTS,
            "Equilibrium constraint requires at least one reactant");
      }
      if (products_.empty())
      {
        throw MicmException(
            MicmSeverity::Error,
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_EMPTY_PRODUCTS,
            "Equilibrium constraint requires at least one product");
      }
      if (equilibrium_constant_ <= 0)
      {
        throw MicmException(
            MicmSeverity::Error,
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_INVALID_EQUILIBRIUM_CONSTANT,
            "Equilibrium constant must be positive");
      }
      for (const auto& r : reactants_)
      {
        if (r.coefficient_ <= 0)
        {
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_INVALID_STOICHIOMETRY,
              "Stoichiometric coefficients must be positive");
        }
      }
      for (const auto& p : products_)
      {
        if (p.coefficient_ <= 0)
        {
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_INVALID_STOICHIOMETRY,
              "Stoichiometric coefficients must be positive");
        }
      }

      // Build species dependencies list (reactants first, then products)
      std::size_t idx = 0;
      for (const auto& r : reactants_)
      {
        species_dependencies_.push_back(r.species_.name_);
        reactant_dependency_indices_.push_back(idx++);
      }
      for (const auto& p : products_)
      {
        species_dependencies_.push_back(p.species_.name_);
        product_dependency_indices_.push_back(idx++);
      }
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    /// @return Species name of the primary algebraic variable
    ///
    /// For equilibrium constraints, we use the first product species as the algebraic row target.
    /// This supports common forms such as K_eq * [B] - [C] = 0 where C is algebraic.
    const std::string& AlgebraicSpecies() const
    {
      return products_[0].species_.name_;
    }

    /// @brief Create function object to compute constraint residual G = K_eq * [reactants]^stoich - [products]^stoich
    ///        Called during solver build (SetConstraintFunctions) to pre-compile residual computation
    /// @param info Constraint information including row index and species indices
    /// @param state_variable_indices Mapping of state variable names to indices
    /// @return Function object that takes (state, forcing) and writes G to forcing[constraint_row]
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)> ResidualFunction(
        const ConstraintInfo& info,
        const auto& state_variable_indices) const
    {
      DenseMatrixPolicy temp_state_variables{ 1, state_variable_indices.size(), 0.0 };

      // Copy data to avoid issues when ConstraintSet is moved
      double K_eq = this->equilibrium_constant_;
      std::vector<double> reactant_stoich;
      std::vector<std::size_t> reactant_state_idx;
      for (std::size_t i = 0; i < this->reactants_.size(); ++i)
      {
        reactant_stoich.push_back(this->reactants_[i].coefficient_);
        reactant_state_idx.push_back(info.state_indices_[reactant_dependency_indices_[i]]);
      }
      std::vector<double> product_stoich;
      std::vector<std::size_t> product_state_idx;
      for (std::size_t i = 0; i < this->products_.size(); ++i)
      {
        product_stoich.push_back(this->products_[i].coefficient_);
        product_state_idx.push_back(info.state_indices_[product_dependency_indices_[i]]);
      }
      std::size_t row_idx = info.row_index_;

      return DenseMatrixPolicy::Function(
          [K_eq, reactant_stoich, reactant_state_idx, product_stoich, product_state_idx, row_idx](auto&& state, auto&& force)
          {
            auto reactant_product = force.GetRowVariable();
            auto product_product = force.GetRowVariable();

            // Initialize reactant_product to K_eq and product_product to 1.0
            state.ForEachRow(
                [K_eq](double& rp, double& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                reactant_product,
                product_product);

            // Multiply in each reactant concentration raised to its stoichiometry
            for (std::size_t i = 0; i < reactant_stoich.size(); ++i)
            {
              const double stoich = reactant_stoich[i];
              const std::size_t species_idx = reactant_state_idx[i];

              state.ForEachRow(
                  [stoich](const double& conc, double& product) { product *= std::pow(std::max(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  reactant_product);
            }

            // Multiply in each product concentration raised to its stoichiometry
            for (std::size_t i = 0; i < product_stoich.size(); ++i)
            {
              const double stoich = product_stoich[i];
              const std::size_t species_idx = product_state_idx[i];

              state.ForEachRow(
                  [stoich](const double& conc, double& product) { product *= std::pow(std::max(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  product_product);
            }

            // Write G = K_eq * [reactants]^stoich - [products]^stoich to the constraint row
            state.ForEachRow(
                [](const double& rp, const double& pp, double& forcing_term) { forcing_term = rp - pp; },
                reactant_product,
                product_product,
                force.GetColumnView(row_idx));
          },
          temp_state_variables,
          temp_state_variables);
    }

    /// @brief Compute Jacobian entries dG/d[species]
    ///        For reactant R with stoichiometry n:
    ///          dG/d[R] = K_eq * n * [R]^(n-1) * prod([other_reactants]^stoich)
    ///        For product P with stoichiometry m:
    ///          dG/d[P] = -m * [P]^(m-1) * prod([other_products]^stoich)
    /// @param info Constraint information including species indices
    /// @param state_variable_indices Mapping of state variable names to indices
    /// @param jacobian_flat_ids Iterator to this constraint's flat Jacobian indices in shared storage
    /// @param jacobian Sparse matrix to store Jacobian values
    /// @return Function object that takes (state_variables, jacobian) and computes partials
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const ConstraintInfo& info,
        const auto& state_variable_indices,
        auto jacobian_flat_ids,
        SparseMatrixPolicy& jacobian) const
    {
      DenseMatrixPolicy temp_state_variables{ 1, state_variable_indices.size(), 0.0 };

      // Pre-compute flat IDs and store them in a vector
      // This avoids iterator issues when the lambda is called multiple times (once per block)
      std::vector<std::size_t> flat_ids;
      flat_ids.reserve(reactants_.size() + products_.size());
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        flat_ids.push_back(*jacobian_flat_ids++);
      }
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        flat_ids.push_back(*jacobian_flat_ids++);
      }

      // Copy data to avoid issues when ConstraintSet is moved
      double K_eq = this->equilibrium_constant_;
      std::vector<double> reactant_stoich;
      std::vector<std::size_t> reactant_state_idx;
      for (std::size_t i = 0; i < this->reactants_.size(); ++i)
      {
        reactant_stoich.push_back(this->reactants_[i].coefficient_);
        reactant_state_idx.push_back(info.state_indices_[reactant_dependency_indices_[i]]);
      }
      std::vector<double> product_stoich;
      std::vector<std::size_t> product_state_idx;
      for (std::size_t i = 0; i < this->products_.size(); ++i)
      {
        product_stoich.push_back(this->products_[i].coefficient_);
        product_state_idx.push_back(info.state_indices_[product_dependency_indices_[i]]);
      }

      return SparseMatrixPolicy::Function(
          [K_eq, reactant_stoich, reactant_state_idx, product_stoich, product_state_idx, flat_ids](
              auto&& state, auto&& jacobian_values)
          {
            // Create temporary variables for computing partials
            auto reactant_product = jacobian_values.GetBlockVariable();
            auto product_product = jacobian_values.GetBlockVariable();
            auto partial_derivative = jacobian_values.GetBlockVariable();

            jacobian_values.ForEachBlock(
                [K_eq](double& rp, double& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                reactant_product,
                product_product);

            for (std::size_t i = 0; i < reactant_stoich.size(); ++i)
            {
              const double stoich = reactant_stoich[i];
              const std::size_t species_idx = reactant_state_idx[i];

              jacobian_values.ForEachBlock(
                  [stoich](const double& conc, double& product) { product *= std::pow(std::max(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  reactant_product);
            }

            for (std::size_t i = 0; i < product_stoich.size(); ++i)
            {
              const double stoich = product_stoich[i];
              const std::size_t species_idx = product_state_idx[i];

              jacobian_values.ForEachBlock(
                  [stoich](const double& conc, double& product) { product *= std::pow(std::max(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  product_product);
            }

            // Compute Jacobian entries for each reactant: dG/d[R_i] = K_eq * stoich_i * prod([other_reactants]^stoich) *
            // [R_i]^(stoich_i-1)
            for (std::size_t i = 0; i < reactant_stoich.size(); ++i)
            {
              const double stoich_i = reactant_stoich[i];
              const std::size_t species_idx_i = reactant_state_idx[i];

              // Compute product of K_eq * all reactants except R_i
              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock([K_eq](double& prod) { prod = K_eq; }, partial_product);

              for (std::size_t j = 0; j < reactant_stoich.size(); ++j)
              {
                if (j != i)  // Skip current species
                {
                  const double stoich_j = reactant_stoich[j];
                  const std::size_t species_idx_j = reactant_state_idx[j];

                  jacobian_values.ForEachBlock(
                      [stoich_j](const double& conc, double& prod) { prod *= std::pow(std::max(0.0, conc), stoich_j); },
                      state.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              // Multiply by stoich_i * [R_i]^(stoich_i-1)
              jacobian_values.ForEachBlock(
                  [stoich_i](const double& conc, const double& prod, double& partial)
                  {
                    if (stoich_i == 1.0)
                    {
                      partial = prod;
                    }
                    else if (conc > 0.0)
                    {
                      partial = stoich_i * prod * std::pow(conc, stoich_i - 1.0);
                    }
                    else
                    {
                      partial = 0.0;
                    }
                  },
                  state.GetConstColumnView(species_idx_i),
                  partial_product,
                  partial_derivative);

              // Subtract partial from Jacobian (matching the SubtractJacobianTerms convention)
              // Use pre-computed flat ID for this reactant
              jacobian_values.ForEachBlock(
                  [](const double& partial, double& jac) { jac -= partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(flat_ids[i]));
            }

            // Compute Jacobian entries for each product: dG/d[P_i] = -stoich_i * prod([other_products]^stoich) *
            // [P_i]^(stoich_i-1)
            for (std::size_t i = 0; i < product_stoich.size(); ++i)
            {
              const double stoich_i = product_stoich[i];
              const std::size_t species_idx_i = product_state_idx[i];

              // Compute product of all products except P_i
              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock([](double& prod) { prod = 1.0; }, partial_product);

              for (std::size_t j = 0; j < product_stoich.size(); ++j)
              {
                if (j != i)  // Skip current species
                {
                  const double stoich_j = product_stoich[j];
                  const std::size_t species_idx_j = product_state_idx[j];

                  jacobian_values.ForEachBlock(
                      [stoich_j](const double& conc, double& prod) { prod *= std::pow(std::max(0.0, conc), stoich_j); },
                      state.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              // Multiply by stoich_i * [P_i]^(stoich_i-1)
              jacobian_values.ForEachBlock(
                  [stoich_i](const double& conc, const double& prod, double& partial)
                  {
                    if (stoich_i == 1.0)
                    {
                      partial = prod;
                    }
                    else if (conc > 0.0)
                    {
                      partial = stoich_i * prod * std::pow(conc, stoich_i - 1.0);
                    }
                    else
                    {
                      partial = 0.0;
                    }
                  },
                  state.GetConstColumnView(species_idx_i),
                  partial_product,
                  partial_derivative);

              // Add partial to Jacobian (note: G = ... - [products], so derivative gets positive sign after subtraction)
              // Use pre-computed flat ID for this product (products come after reactants in flat_ids)
              jacobian_values.ForEachBlock(
                  [](const double& partial, double& jac) { jac += partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(flat_ids[reactant_stoich.size() + i]));
            }
          },
          temp_state_variables,
          jacobian);
    }
  };

}  // namespace micm
