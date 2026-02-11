// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/equilibrium_constraint.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>
#include <concepts>
#include <type_traits>

namespace micm
{

  /// @brief This class uses std::variant to hold different constraint types.
  ///        Each constraint provides:
  ///        - A residual function G(y) that should equal zero when the constraint is satisfied
  ///        - Jacobian entries dG/dy for each species the constraint depends on
  class Constraint
  {
   public:
    using ConstraintVariant = std::variant<EquilibriumConstraint>;

    ConstraintVariant constraint_;

    template<typename T>
      requires std::same_as<std::decay_t<T>, EquilibriumConstraint>
    Constraint(T&& constraint)
        : constraint_(std::forward<T>(constraint))
    {
    }

    /// @brief Get the constraint name
    /// @return Constraint name
    std::string GetName() const
    {
      return std::visit([](const auto& c) { return c.name_; }, constraint_);
    }

    /// @brief Get species dependencies
    /// @return Vector of species names this constraint depends on
    const std::vector<std::string>& GetSpeciesDependencies() const
    {
      return std::visit([](const auto& c) -> const std::vector<std::string>& { return c.species_dependencies_; }, constraint_);
    }

    /// @brief Get the number of species this constraint depends on
    /// @return Number of dependent species
    std::size_t NumberOfDependencies() const
    {
      return std::visit([](const auto& c) { return c.species_dependencies_.size(); }, constraint_);
    }

    /// @brief Evaluate the constraint residual G(y)
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to positions in concentrations
    /// @return Residual value (should be 0 when constraint is satisfied)
    double Residual(const double* concentrations, const std::size_t* indices) const
    {
      return std::visit([&](const auto& c) { return c.Residual(concentrations, indices); }, constraint_);
    }

    /// @brief Compute partial derivatives dG/d[species] for each dependent species
    /// @param concentrations Pointer to species concentrations (row of state matrix)
    /// @param indices Pointer to indices mapping species_dependencies_ to positions in concentrations
    /// @param jacobian Output buffer for partial derivatives dG/d[species] (same order as species_dependencies_)
    void Jacobian(const double* concentrations, const std::size_t* indices, double* jacobian) const
    {
      std::visit([&](const auto& c) { c.Jacobian(concentrations, indices, jacobian); }, constraint_);
    }

    /// @brief Returns the species whose state row should be replaced by this algebraic constraint
    /// @return Algebraic species name
    const std::string& GetAlgebraicSpecies() const
    {
      return std::visit([](const auto& c) -> const std::string& { return c.AlgebraicSpecies(); }, constraint_);
    }
  };

}  // namespace micm
