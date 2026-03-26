// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <concepts>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>
#include <unordered_set>

namespace micm
{

  /// @brief This class uses std::variant to hold different constraint types.
  ///        Each constraint provides:
  ///        - A residual function G(y) that should equal zero when the constraint is satisfied
  ///        - Jacobian entries dG/dy for each species the constraint depends on
  class Constraint
  {
   public:
    using ConstraintVariant = std::variant<EquilibriumConstraint, LinearConstraint>;

    ConstraintVariant constraint_;

    template<typename T>
      requires std::same_as<std::decay_t<T>, EquilibriumConstraint> || std::same_as<std::decay_t<T>, LinearConstraint>
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

    /// @brief Get the custom paramter names
    /// @return A set of parameter names
    std::unordered_set<std::string> GetParameterNames() const
    {
      return std::visit([](const auto& c) { return c.parameters_; }, constraint_);
    }

    /// @brief Returns the species whose state row should be replaced by this algebraic constraint
    /// @return Algebraic species name
    const std::string& AlgebraicSpecies() const
    {
      return std::visit([](const auto& c) -> const std::string& { return c.AlgebraicSpecies(); }, constraint_);
    }

    /// @brief Get species dependencies
    /// @return Vector of species names this constraint depends on
    const std::vector<std::string>& SpeciesDependencies() const
    {
      return std::visit(
          [](const auto& c) -> const std::vector<std::string>& { return c.species_dependencies_; }, constraint_);
    }

    /// @brief Get the number of species this constraint depends on
    /// @return Number of dependent species
    std::size_t NumberOfDependencies() const
    {
      return std::visit([](const auto& c) { return c.species_dependencies_.size(); }, constraint_);
    }

    /// @brief Get a function object to compute the constraint residual
    ///        This returns a reusable function that can be invoked multiple times
    /// @param info Constraint information including species indices and row index
    /// @param state_variable_indices Map from species names to state variable indices
    /// @return Function object that takes (state_variables, forcing) and computes the residual
    template<typename DenseMatrixPolicy>
    auto ResidualFunction(const ConstraintInfo& info, const auto& state_variable_indices) const
    {
      return std::visit(
          [&info, &state_variable_indices](const auto& c)
          { return c.template ResidualFunction<DenseMatrixPolicy>(info, state_variable_indices); },
          constraint_);
    }

    /// @brief Get a function object to compute the constraint Jacobian
    ///        This returns a reusable function that can be invoked multiple times
    /// @param info Constraint information including species indices and Jacobian flat IDs
    /// @param state_variable_indices Map from species names to state variable indices
    /// @param jacobian_flat_ids Iterator to the jacobian flat IDs for this constraint
    /// @param jacobian Sparse matrix to store Jacobian values
    /// @return Function object that takes (state_variables, jacobian) and computes partials
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    auto JacobianFunction(
        const ConstraintInfo& info,
        const auto& state_variable_indices,
        auto jacobian_flat_ids,
        SparseMatrixPolicy& jacobian) const
    {
      return std::visit(
          [&info, &state_variable_indices, jacobian_flat_ids, &jacobian](const auto& c)
          {
            return c.template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
                info, state_variable_indices, jacobian_flat_ids, jacobian);
          },
          constraint_);
    }
  };

}  // namespace micm
