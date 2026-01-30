// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace micm
{

  /// @brief Abstract base class for algebraic constraints G(y) = 0
  ///
  /// Constraints define algebraic relations that must be satisfied by the species
  /// concentrations. They are used in DAE (Differential-Algebraic Equation) solvers
  /// to enforce conditions like chemical equilibrium or mass conservation.
  ///
  /// Each constraint provides:
  /// - A residual function G(y) that should equal zero when the constraint is satisfied
  /// - Jacobian entries dG/dy for each species the constraint depends on
  class Constraint
  {
   public:
    /// @brief Name of the constraint (for identification/debugging)
    std::string name_;

    /// @brief Names of species this constraint depends on
    std::vector<std::string> species_dependencies_;

    /// @brief Default constructor
    Constraint() = default;

    /// @brief Constructor with name
    /// @param name Constraint identifier
    explicit Constraint(const std::string& name)
        : name_(name)
    {
    }

    /// @brief Virtual destructor
    virtual ~Constraint() = default;

    /// @brief Evaluate the constraint residual G(y)
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to positions in concentrations
    /// @return Residual value (should be 0 when constraint is satisfied)
    virtual double Residual(const double* concentrations, const std::size_t* indices) const = 0;

    /// @brief Compute partial derivatives dG/d[species] for each dependent species
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to positions in concentrations
    /// @param jacobian Output buffer for partial derivatives dG/d[species] (same order as species_dependencies_)
    virtual void Jacobian(const double* concentrations, const std::size_t* indices, double* jacobian) const = 0;

    /// @brief Get the number of species this constraint depends on
    /// @return Number of dependent species
    std::size_t NumberOfDependencies() const
    {
      return species_dependencies_.size();
    }
  };

}  // namespace micm
