// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/util/types.hpp>

#include <concepts>
#include <cstddef>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  /// @brief This class uses std::variant to hold different constraint types.
  ///        Each constraint provides:
  ///        - A residual function G(y) that should equal zero when the constraint is satisfied
  ///        - Jacobian entries dG/dy for each species the constraint depends on
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class Constraint
  {
   public:
    using ConstraintVariant = std::variant<
        EquilibriumConstraint<DenseMatrixPolicy, SparseMatrixPolicy>,
        LinearConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>;

    ConstraintVariant constraint_;

    template<typename T>
      requires std::same_as<std::decay_t<T>, EquilibriumConstraint<DenseMatrixPolicy, SparseMatrixPolicy>> ||
               std::same_as<std::decay_t<T>, LinearConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>
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

    /// @brief Get the custom parameter names
    /// @return A set of parameter names
    std::vector<std::string> GetParameterNames() const
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
    Index NumberOfDependencies() const
    {
      return std::visit([](const auto& c) { return c.species_dependencies_.size(); }, constraint_);
    }

    /// @brief Apply constraint parameter update for all grid cells (e.g., temperature-dependent K_eq)
    ///        Called directly from ConstraintSet::UpdateStateParameters.
    void ApplyConstraintParameter(
        const ConstraintInfo& info,
        const typename DenseMatrixPolicy::template VectorType<Conditions>& conditions,
        DenseMatrixPolicy& state_param) const
    {
      if (auto* eq = std::get_if<EquilibriumConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        eq->ApplyConstraintParameter(info, conditions, state_param);
      }
      else if (auto* lin = std::get_if<LinearConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        lin->ApplyConstraintParameter(info, conditions, state_param);
      }
    }

    void SetStateIndices(const ConstraintInfo& info, auto& jacobian_flat_ids)
    {
      return std::visit([&info, &jacobian_flat_ids](auto& c) { c.SetStateIndices(info, jacobian_flat_ids); }, constraint_);
    }

    /// @brief Add constraint residual G to forcing vector for all grid cells
    ///        Called directly from ConstraintSet::AddForcingTerms.
    void AddResidual(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        DenseMatrixPolicy& forcing) const
    {
      if (auto* eq = std::get_if<EquilibriumConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        eq->AddResidual(info, state, state_param, forcing);
      }
      else if (auto* lin = std::get_if<LinearConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        lin->AddResidual(info, state, state_param, forcing);
      }
    }

    /// @brief Subtract Jacobian partial derivatives from Jacobian matrix for all grid cells
    ///        Called directly from ConstraintSet::SubtractJacobianTerms.
    void SubtractJacobian(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        SparseMatrixPolicy& jacobian) const
    {
      if (auto* eq = std::get_if<EquilibriumConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        eq->SubtractJacobian(info, state, state_param, jacobian);
      }
      else if (auto* lin = std::get_if<LinearConstraint<DenseMatrixPolicy, SparseMatrixPolicy>>(&constraint_))
      {
        lin->SubtractJacobian(info, state, state_param, jacobian);
      }
    }
  };

}  // namespace micm
