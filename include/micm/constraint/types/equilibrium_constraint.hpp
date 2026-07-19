// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <system_error>
#include <vector>

namespace micm
{
  /// @brief Define parameters for Van't Hoff equation
  struct VantHoffParam
  {
    Real K_HLC_ref_;                    // Henry’s Law constant at the reference temperature (typically 298.15K)
    Real delta_H_;                      // Enthalpy of dissolution (J/mol)
    Real R_ = constants::GAS_CONSTANT;  // (J/mol·K)
    Real T_ref_ = 298.15;
  };

  /// @brief Constraint for chemical equilibrium with temperature-dependent K_eq using Van't Hoff equation
  ///        For a reversible reaction: aA + bB <-> cC + dD
  ///        The equilibrium constraint is: G = K_eq(T) * [A]^a * [B]^b - [C]^c * [D]^d = 0
  ///        where K_eq(T) = K_HLC_ref * exp((delta_H / R) * (1/T - 1/T_ref))
  class EquilibriumConstraint
  {
   public:
    /// @brief Name of the constraint, used when generating state parameter name
    std::string name_;

    /// @brief Algebraic species
    Species algebraic_species_;

    /// @brief Names of species this constraint depends on
    std::vector<std::string> species_dependencies_;

    /// @brief Reactant species and their stoichiometric coefficients
    std::vector<StoichSpecies> reactants_;

    /// @brief Product species and their stoichiometric coefficients
    std::vector<StoichSpecies> products_;

    /// @brief For equilibrium constraints, this contains a single parameter K_eq
    std::vector<std::string> parameters_;

    /// @brief Temperature-dependent Henry’s Law Constant
    std::function<Real(const Conditions&)> equilibrium_constant_function_;

   private:
    /// @brief Van't Hoff equation parameter used to calculate Henry’s Law Constant
    VantHoffParam vant_hoff_param_;

    /// @brief Indices into the reactants_ vector for each species dependency
    std::vector<Index> reactant_dependency_indices_;

    /// @brief Indices into the products_ vector for each species dependency
    std::vector<Index> product_dependency_indices_;

   public:
    /// @brief Default constructor
    EquilibriumConstraint() = default;

    /// @brief Construct an equilibrium constraint.
    ///        Validates that equilibrium constraint > 0.
    ///        Builds species_dependencies_ by concatenating reactants then products.
    ///        Stores index mappings for efficient Jacobian computation.
    ///        Stores a temperature-dependent equilibrium constant function.
    /// @param name Constraint identifier
    /// @param algebraic_species Species whose row is replaced by this algebraic constraint
    /// @param reactants Vector of StoichSpecies (species, stoichiometry) for reactants
    /// @param products Vector of StoichSpecies (species, stoichiometry) for products
    /// @param vant_hoff_param Parameters for Van't Hoff equation
    EquilibriumConstraint(
        const std::string& name,
        const Species& algebraic_species,
        std::vector<StoichSpecies> reactants,
        std::vector<StoichSpecies> products,
        VantHoffParam vant_hoff_param)
        : name_(name),
          algebraic_species_(algebraic_species),
          reactants_(std::move(reactants)),
          products_(std::move(products)),
          vant_hoff_param_(vant_hoff_param)
    {
      if (reactants_.empty())
      {
        throw MicmException(
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_EMPTY_REACTANTS,
            "Equilibrium constraint requires at least one reactant");
      }
      if (products_.empty())
      {
        throw MicmException(
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_EMPTY_PRODUCTS,
            "Equilibrium constraint requires at least one product");
      }
      for (const auto& r : reactants_)
      {
        if (r.coefficient_ <= 0)
        {
          throw MicmException(
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
              MICM_ERROR_CATEGORY_CONSTRAINT,
              MICM_CONSTRAINT_ERROR_CODE_INVALID_STOICHIOMETRY,
              "Stoichiometric coefficients must be positive");
        }
      }
      if (vant_hoff_param_.K_HLC_ref_ <= 0)
      {
        throw MicmException(
            MICM_ERROR_CATEGORY_CONSTRAINT,
            MICM_CONSTRAINT_ERROR_CODE_INVALID_EQUILIBRIUM_CONSTANT,
            "Henry’s Law constant (K_HLC_ref) must be positive");
      }
      parameters_.push_back(name);

      // Build species dependencies list (reactants first, then products)
      Index idx = 0;
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

      equilibrium_constant_function_ = [p = vant_hoff_param_](const Conditions& condition)
      {
        Real T = condition.temperature_;
        return p.K_HLC_ref_ * std::exp((p.delta_H_ / p.R_) * (1.0 / T - 1.0 / p.T_ref_));
      };
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    /// @return Species name of the explicitly set algebraic variable
    const std::string& AlgebraicSpecies() const
    {
      return algebraic_species_.name_;
    }

    /// @brief Create function object to update temperature-dependent K_eq parameter
    ///        Returns a function that computes K_eq(T) for each grid cell using Van't Hoff equation
    ///        Called during solver's UpdateStateParameters phase before each solve
    /// @param info Constraint information including state parameter indices
    /// @return Function object that takes (conditions, state_param) and writes K_eq(T) to state_param[K_eq_idx]
    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)> ConstraintParameterFunction(
        const ConstraintInfo& info) const
    {
      Index K_eq_idx = info.state_param_indices_[0];  // equilibrium constant index

      return [K_eq_idx, eq_func = equilibrium_constant_function_](
                 const std::vector<Conditions>& conditions, DenseMatrixPolicy& state_param)
      {
        // For each grid cell, compute K_eq at current temperature
        state_param.ForEachRow(
            [eq_func](const Conditions& cond, Real& K_eq) { K_eq = eq_func(cond); },
            conditions,
            state_param.GetColumnView(K_eq_idx));
      };
    }

    /// @brief Create function object to compute equilibrium constraint residual for all grid cells
    ///        Computes G = K_eq(T) * prod([reactants]^stoich) - prod([products]^stoich) for the algebraic constraint
    ///        Called during solver build (SetConstraintFunctions) to pre-compile residual computation
    /// @param info Constraint information including row index, species indices, and parameter indices
    /// @param state_variable_indices Mapping of state variable names to column indices in state matrix
    /// @param state_parameter_indices Mapping of state parameter names to column indices in state_param matrix
    /// @return Function object that takes (state, state_param, forcing) and writes residual G to forcing[constraint_row]
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ResidualFunction(
        const ConstraintInfo& info,
        const auto& state_variable_indices,
        const auto& state_parameter_indices) const
    {
      // Copy data to avoid issues when ConstraintSet is moved
      std::vector<Real> reactant_stoich;
      std::vector<Index> reactant_state_idx;
      for (Index i = 0; i < this->reactants_.size(); ++i)
      {
        reactant_stoich.push_back(this->reactants_[i].coefficient_);
        reactant_state_idx.push_back(info.state_indices_[reactant_dependency_indices_[i]]);
      }
      std::vector<Real> product_stoich;
      std::vector<Index> product_state_idx;
      for (Index i = 0; i < this->products_.size(); ++i)
      {
        product_stoich.push_back(this->products_[i].coefficient_);
        product_state_idx.push_back(info.state_indices_[product_dependency_indices_[i]]);
      }
      Index row_idx = info.row_index_;
      Index K_eq_idx = info.state_param_indices_[0];  // contains only one parameter (equilibrium constant)

      DenseMatrixPolicy temp_state_variables{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy temp_state_parameters{ 1, state_parameter_indices.size(), 0.0 };

      return DenseMatrixPolicy::Function(
          [=](auto&& state, auto&& state_param, auto&& force)
          {
            auto reactant_product = force.GetRowVariable();
            auto product_product = force.GetRowVariable();

            // Initialize reactant_product to K_eq and product_product to 1.0
            state.ForEachRow(
                [](const Real& K_eq, Real& rp, Real& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                state_param.GetConstColumnView(K_eq_idx),
                reactant_product,
                product_product);

            // Multiply in each reactant concentration raised to its stoichiometry
            for (Index i = 0; i < reactant_stoich.size(); ++i)
            {
              const Real stoich = reactant_stoich[i];
              const Index species_idx = reactant_state_idx[i];

              state.ForEachRow(
                  [stoich](const Real& conc, Real& product) { product *= std::pow(std::max<Real>(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  reactant_product);
            }

            // Multiply in each product concentration raised to its stoichiometry
            for (Index i = 0; i < product_stoich.size(); ++i)
            {
              const Real stoich = product_stoich[i];
              const Index species_idx = product_state_idx[i];

              state.ForEachRow(
                  [stoich](const Real& conc, Real& product) { product *= std::pow(std::max<Real>(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  product_product);
            }

            // Write G = K_eq * [reactants]^stoich - [products]^stoich to the constraint row
            state.ForEachRow(
                [](const Real& rp, const Real& pp, Real& forcing_term) { forcing_term = rp - pp; },
                reactant_product,
                product_product,
                force.GetColumnView(row_idx));
          },
          temp_state_variables,
          temp_state_parameters,
          temp_state_variables);
    }

    /// @brief Create function object to compute Jacobian partial derivatives dG/d[species] for all grid cells
    ///        For reactant R with stoichiometry n:
    ///          dG/d[R] = K_eq(T) * n * [R]^(n-1) * prod([other_reactants]^stoich)
    ///        For product P with stoichiometry m:
    ///          dG/d[P] = -m * [P]^(m-1) * prod([other_products]^stoich)
    ///        Called during solver build (SetConstraintFunctions) to pre-compile Jacobian computation
    /// @param info Constraint information including row index, species indices, and parameter indices
    /// @param state_variable_indices Mapping of state variable names to column indices in state matrix
    /// @param state_parameter_indices Mapping of state parameter names to column indices in state_param matrix
    /// @param jacobian_flat_ids Iterator to this constraint's flat Jacobian indices in sparse matrix storage
    /// @param jacobian Sparse matrix reference (used for type information)
    /// @return Function object that takes (state, state_param, jacobian_values) and writes partials to sparse Jacobian
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const ConstraintInfo& info,
        const auto& state_variable_indices,
        const auto& state_parameter_indices,
        auto jacobian_flat_ids,
        SparseMatrixPolicy& jacobian) const
    {
      // Pre-compute flat IDs and store them in a vector
      // This avoids iterator issues when the lambda is called multiple times (once per block)
      std::vector<Index> flat_ids;
      flat_ids.reserve(reactants_.size() + products_.size());
      for (Index i = 0; i < reactants_.size(); ++i)
      {
        flat_ids.push_back(*jacobian_flat_ids++);
      }
      for (Index i = 0; i < products_.size(); ++i)
      {
        flat_ids.push_back(*jacobian_flat_ids++);
      }

      std::vector<Real> reactant_stoich;
      std::vector<Index> reactant_state_idx;
      for (Index i = 0; i < this->reactants_.size(); ++i)
      {
        reactant_stoich.push_back(this->reactants_[i].coefficient_);
        reactant_state_idx.push_back(info.state_indices_[reactant_dependency_indices_[i]]);
      }
      std::vector<Real> product_stoich;
      std::vector<Index> product_state_idx;
      for (Index i = 0; i < this->products_.size(); ++i)
      {
        product_stoich.push_back(this->products_[i].coefficient_);
        product_state_idx.push_back(info.state_indices_[product_dependency_indices_[i]]);
      }

      Index K_eq_idx = info.state_param_indices_[0];  // contains only one parameter (equilibrium constant)

      DenseMatrixPolicy temp_state_variables{ 1, state_variable_indices.size(), 0.0 };
      DenseMatrixPolicy temp_state_parameters{ 1, state_parameter_indices.size(), 0.0 };

      return SparseMatrixPolicy::Function(
          [=](auto&& state, auto&& state_param, auto&& jacobian_values)
          {
            // Create temporary variables for computing partials
            auto reactant_product = jacobian_values.GetBlockVariable();
            auto product_product = jacobian_values.GetBlockVariable();
            auto partial_derivative = jacobian_values.GetBlockVariable();

            jacobian_values.ForEachBlock(
                [](const Real& K_eq, Real& rp, Real& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                state_param.GetConstColumnView(K_eq_idx),
                reactant_product,
                product_product);

            for (Index i = 0; i < reactant_stoich.size(); ++i)
            {
              const Real stoich = reactant_stoich[i];
              const Index species_idx = reactant_state_idx[i];

              jacobian_values.ForEachBlock(
                  [stoich](const Real& conc, Real& product) { product *= std::pow(std::max<Real>(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  reactant_product);
            }

            for (Index i = 0; i < product_stoich.size(); ++i)
            {
              const Real stoich = product_stoich[i];
              const Index species_idx = product_state_idx[i];

              jacobian_values.ForEachBlock(
                  [stoich](const Real& conc, Real& product) { product *= std::pow(std::max<Real>(0.0, conc), stoich); },
                  state.GetConstColumnView(species_idx),
                  product_product);
            }

            // Compute Jacobian entries for each reactant: dG/d[R_i] = K_eq * stoich_i * prod([other_reactants]^stoich) *
            // [R_i]^(stoich_i-1)
            for (Index i = 0; i < reactant_stoich.size(); ++i)
            {
              const Real stoich_i = reactant_stoich[i];
              const Index species_idx_i = reactant_state_idx[i];

              // Compute product of K_eq * all reactants except R_i
              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock(
                  [](const Real& K_eq, Real& prod) { prod = K_eq; },
                  state_param.GetConstColumnView(K_eq_idx),
                  partial_product);

              for (Index j = 0; j < reactant_stoich.size(); ++j)
              {
                if (j != i)  // Skip current species
                {
                  const Real stoich_j = reactant_stoich[j];
                  const Index species_idx_j = reactant_state_idx[j];

                  jacobian_values.ForEachBlock(
                      [stoich_j](const Real& conc, Real& prod) { prod *= std::pow(std::max<Real>(0.0, conc), stoich_j); },
                      state.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              // Multiply by stoich_i * [R_i]^(stoich_i-1)
              jacobian_values.ForEachBlock(
                  [stoich_i](const Real& conc, const Real& prod, Real& partial)
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
                  [](const Real& partial, Real& jac) { jac -= partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(flat_ids[i]));
            }

            // Compute Jacobian entries for each product: dG/d[P_i] = -stoich_i * prod([other_products]^stoich) *
            // [P_i]^(stoich_i-1)
            for (Index i = 0; i < product_stoich.size(); ++i)
            {
              const Real stoich_i = product_stoich[i];
              const Index species_idx_i = product_state_idx[i];

              // Compute product of all products except P_i
              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock([](Real& prod) { prod = 1.0; }, partial_product);

              for (Index j = 0; j < product_stoich.size(); ++j)
              {
                if (j != i)  // Skip current species
                {
                  const Real stoich_j = product_stoich[j];
                  const Index species_idx_j = product_state_idx[j];

                  jacobian_values.ForEachBlock(
                      [stoich_j](const Real& conc, Real& prod) { prod *= std::pow(std::max<Real>(0.0, conc), stoich_j); },
                      state.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              // Multiply by stoich_i * [P_i]^(stoich_i-1)
              jacobian_values.ForEachBlock(
                  [stoich_i](const Real& conc, const Real& prod, Real& partial)
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
                  [](const Real& partial, Real& jac) { jac += partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(flat_ids[reactant_stoich.size() + i]));
            }
          },
          temp_state_variables,
          temp_state_parameters,
          jacobian);
    }
  };

}  // namespace micm
