// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_error.hpp>
#include <micm/system/stoich_species.hpp>

#include <cmath>
#include <cstddef>
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
    ///        Validates that equilibrium constraint >= 0
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
        throw std::system_error(make_error_code(MicmConstraintErrc::EmptyReactants));
      }
      if (products_.empty())
      {
        throw std::system_error(make_error_code(MicmConstraintErrc::EmptyProducts));
      }
      if (equilibrium_constant_ <= 0)
      {
        throw std::system_error(make_error_code(MicmConstraintErrc::InvalidEquilibriumConstant));
      }
      for (const auto& r : reactants_)
      {
        if (r.coefficient_ <= 0)
        {
          throw std::system_error(make_error_code(MicmConstraintErrc::InvalidStoichiometry));
        }
      }
      for (const auto& p : products_)
      {
        if (p.coefficient_ <= 0)
        {
          throw std::system_error(make_error_code(MicmConstraintErrc::InvalidStoichiometry));
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

    /// @brief Evaluate the equilibrium constraint residual
    ///        G = K_eq * prod([reactants]^stoich) - prod([products]^stoich)
    ///        At equilibrium, G = 0
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to concentrations
    /// @return Residual value
    double Residual(const double* concentrations, const std::size_t* indices) const
    {
      // Compute product of reactant concentrations raised to stoichiometric powers
      double reactant_product = 1.0;
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        reactant_product *= std::pow(conc, reactants_[i].coefficient_);
      }

      // Compute product of product concentrations raised to stoichiometric powers
      double product_product = 1.0;
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        double conc = concentrations[indices[product_dependency_indices_[i]]];
        product_product *= std::pow(conc, products_[i].coefficient_);
      }

      // G = K_eq * [reactants] - [products]
      return equilibrium_constant_ * reactant_product - product_product;
    }

    /// @brief Compute Jacobian entries dG/d[species]
    ///        For reactant R with stoichiometry n:
    ///          dG/d[R] = K_eq * n * [R]^(n-1) * prod([other_reactants]^stoich)
    ///        For product P with stoichiometry m:
    ///          dG/d[P] = -m * [P]^(m-1) * prod([other_products]^stoich)
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to concentrations
    /// @param jacobian Output buffer for partial derivatives (same order as species_dependencies_)
    void Jacobian(const double* concentrations, const std::size_t* indices, double* jacobian) const
    {
      // Initialize jacobian entries to zero
      for (std::size_t i = 0; i < species_dependencies_.size(); ++i)
      {
        jacobian[i] = 0.0;
      }

      // Compute full reactant and product terms
      double reactant_product = 1.0;
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        reactant_product *= std::pow(conc, reactants_[i].coefficient_);
      }

      double product_product = 1.0;
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        double conc = concentrations[indices[product_dependency_indices_[i]]];
        product_product *= std::pow(conc, products_[i].coefficient_);
      }

      // Jacobian for reactants: dG/d[R] = K_eq * n * [R]^(n-1) * prod([others])
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        double stoich = reactants_[i].coefficient_;

        if (conc > 0)
        {
          // dG/d[R] = K_eq * stoich * reactant_product / [R]
          jacobian[reactant_dependency_indices_[i]] = equilibrium_constant_ * stoich * reactant_product / conc;
        }
        else if (stoich == 1.0)
        {
          // Special case: if conc = 0 and stoich = 1, derivative is K_eq * prod(others)
          // Recompute product of other reactants, scaled by K_eq
          double others = equilibrium_constant_;
          for (std::size_t j = 0; j < reactants_.size(); ++j)
          {
            if (j != i)
            {
              double other_conc = concentrations[indices[reactant_dependency_indices_[j]]];
              others *= std::pow(other_conc, reactants_[j].coefficient_);
            }
          }
          jacobian[reactant_dependency_indices_[i]] = others;
        }
        // else: derivative is 0 when conc = 0 and stoich > 1
      }

      // Jacobian for products: dG/d[P] = -m * [P]^(m-1) * prod([others])
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        double conc = concentrations[indices[product_dependency_indices_[i]]];
        double stoich = products_[i].coefficient_;

        if (conc > 0)
        {
          // dG/d[P] = -stoich * product_product / [P]
          jacobian[product_dependency_indices_[i]] = -stoich * product_product / conc;
        }
        else if (stoich == 1.0)
        {
          // Special case: if conc = 0 and stoich = 1, derivative is -prod(others)
          double others = -1.0;
          for (std::size_t j = 0; j < products_.size(); ++j)
          {
            if (j != i)
            {
              double other_conc = concentrations[indices[product_dependency_indices_[j]]];
              others *= std::pow(other_conc, products_[j].coefficient_);
            }
          }
          jacobian[product_dependency_indices_[i]] = others;
        }
        // else: derivative is 0 when conc = 0 and stoich > 1
      }
    }
  };

}  // namespace micm
