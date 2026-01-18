// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Constraint for chemical equilibrium: K_eq = [products]^stoich / [reactants]^stoich
  ///
  /// For a reversible reaction: aA + bB <-> cC + dD
  /// The equilibrium constraint is: G = K_eq * [A]^a * [B]^b - [C]^c * [D]^d = 0
  ///
  /// This can also be written in terms of forward/backward rate constants:
  /// G = k_f * [A]^a * [B]^b - k_b * [C]^c * [D]^d = 0
  /// where K_eq = k_f / k_b
  class EquilibriumConstraint : public Constraint
  {
   public:
    /// @brief Reactant species names and their stoichiometric coefficients
    std::vector<std::pair<std::string, double>> reactants_;

    /// @brief Product species names and their stoichiometric coefficients
    std::vector<std::pair<std::string, double>> products_;

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
    /// @param name Constraint identifier
    /// @param reactants Vector of (species_name, stoichiometry) pairs for reactants
    /// @param products Vector of (species_name, stoichiometry) pairs for products
    /// @param equilibrium_constant K_eq = [products]/[reactants] at equilibrium
    EquilibriumConstraint(
        const std::string& name,
        const std::vector<std::pair<std::string, double>>& reactants,
        const std::vector<std::pair<std::string, double>>& products,
        double equilibrium_constant)
        : Constraint(name),
          reactants_(reactants),
          products_(products),
          equilibrium_constant_(equilibrium_constant)
    {
      if (equilibrium_constant_ <= 0)
      {
        throw std::invalid_argument("Equilibrium constant must be positive");
      }

      // Build species dependencies list (reactants first, then products)
      std::size_t idx = 0;
      for (const auto& r : reactants_)
      {
        species_dependencies_.push_back(r.first);
        reactant_dependency_indices_.push_back(idx++);
      }
      for (const auto& p : products_)
      {
        species_dependencies_.push_back(p.first);
        product_dependency_indices_.push_back(idx++);
      }
    }

    /// @brief Evaluate the equilibrium constraint residual
    ///
    /// G = K_eq * prod([reactants]^stoich) - prod([products]^stoich)
    ///
    /// At equilibrium, G = 0
    ///
    /// @param concentrations Vector of all species concentrations
    /// @param indices Indices mapping species_dependencies_ to concentrations vector
    /// @return Residual value
    double Residual(const std::vector<double>& concentrations,
                    const std::vector<std::size_t>& indices) const override
    {
      // Compute product of reactant concentrations raised to stoichiometric powers
      double reactant_product = 1.0;
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        reactant_product *= std::pow(conc, reactants_[i].second);
      }

      // Compute product of product concentrations raised to stoichiometric powers
      double product_product = 1.0;
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        double conc = concentrations[indices[product_dependency_indices_[i]]];
        product_product *= std::pow(conc, products_[i].second);
      }

      // G = K_eq * [reactants] - [products]
      return equilibrium_constant_ * reactant_product - product_product;
    }

    /// @brief Compute Jacobian entries dG/d[species]
    ///
    /// For reactant R with stoichiometry n:
    ///   dG/d[R] = K_eq * n * [R]^(n-1) * prod([other_reactants]^stoich)
    ///
    /// For product P with stoichiometry m:
    ///   dG/d[P] = -m * [P]^(m-1) * prod([other_products]^stoich)
    ///
    /// @param concentrations Vector of all species concentrations
    /// @param indices Indices mapping species_dependencies_ to concentrations vector
    /// @return Vector of partial derivatives in same order as species_dependencies_
    std::vector<double> Jacobian(const std::vector<double>& concentrations,
                                 const std::vector<std::size_t>& indices) const override
    {
      std::vector<double> jacobian(species_dependencies_.size(), 0.0);

      // Compute full reactant and product terms
      double reactant_product = 1.0;
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        reactant_product *= std::pow(conc, reactants_[i].second);
      }

      double product_product = 1.0;
      for (std::size_t i = 0; i < products_.size(); ++i)
      {
        double conc = concentrations[indices[product_dependency_indices_[i]]];
        product_product *= std::pow(conc, products_[i].second);
      }

      // Jacobian for reactants: dG/d[R] = K_eq * n * [R]^(n-1) * prod([others])
      for (std::size_t i = 0; i < reactants_.size(); ++i)
      {
        double conc = concentrations[indices[reactant_dependency_indices_[i]]];
        double stoich = reactants_[i].second;

        if (conc > 0)
        {
          // dG/d[R] = K_eq * stoich * reactant_product / [R]
          jacobian[reactant_dependency_indices_[i]] = equilibrium_constant_ * stoich * reactant_product / conc;
        }
        else if (stoich == 1.0)
        {
          // Special case: if conc = 0 and stoich = 1, derivative is K_eq * prod(others)
          double others = reactant_product;  // This is 0 if conc = 0, need to recompute
          others = equilibrium_constant_;
          for (std::size_t j = 0; j < reactants_.size(); ++j)
          {
            if (j != i)
            {
              double other_conc = concentrations[indices[reactant_dependency_indices_[j]]];
              others *= std::pow(other_conc, reactants_[j].second);
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
        double stoich = products_[i].second;

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
              others *= std::pow(other_conc, products_[j].second);
            }
          }
          jacobian[product_dependency_indices_[i]] = others;
        }
        // else: derivative is 0 when conc = 0 and stoich > 1
      }

      return jacobian;
    }
  };

}  // namespace micm
