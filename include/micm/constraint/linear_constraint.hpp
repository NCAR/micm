// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_error.hpp>
#include <micm/system/stoich_species.hpp>

#include <cstddef>
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
      for (const auto& term : terms_)
      {
        species_dependencies_.push_back(term.species_.name_);
      }
    }

    /// @brief Evaluate the linear constraint residual
    ///        G = sum(coeff[i] * [species[i]]) - constant
    ///        When satisfied, G = 0
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to concentrations
    /// @return Residual value
    double Residual(const double* concentrations, const std::size_t* indices) const
    {
      double sum = 0.0;
      for (std::size_t i = 0; i < terms_.size(); ++i)
      {
        sum += terms_[i].coefficient_ * concentrations[indices[i]];
      }
      return sum - constant_;
    }

    /// @brief Compute Jacobian entries dG/d[species]
    ///        For a linear constraint, the Jacobian is simply the coefficients:
    ///        dG/d[species[i]] = coeff[i]
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to concentrations
    /// @param jacobian Output buffer for partial derivatives (same order as species_dependencies_)
    void Jacobian(const double* concentrations, const std::size_t* indices, double* jacobian) const
    {
      for (std::size_t i = 0; i < terms_.size(); ++i)
      {
        jacobian[i] = terms_[i].coefficient_;
      }
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    ///        For linear constraints, the last species is used in the terms list as the algebraic row target.
    /// @return Species name of the primary algebraic variable
    const std::string& AlgebraicSpecies() const
    {
      return terms_.back().species_.name_;
    }
  };

}  // namespace micm
