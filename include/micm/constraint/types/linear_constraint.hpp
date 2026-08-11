// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/types.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Constraint for linear relationships: sum(coeff[i] * [species[i]]) = constant
  ///        For example: A + B + C = 1.0 represents a conservation law
  ///        The linear constraint is: G = c1*[A] + c2*[B] + c3*[C] - constant = 0
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class LinearConstraint
  {
    using DenseMatrix = DenseMatrixPolicy;
    using SparseMatrix = SparseMatrixPolicy;
    template<class U>
    using Vector = typename DenseMatrix::template VectorType<U>;
    template<class U>
    using VectorView = typename Vector<U>::ConstViewType;
    template<class U>
    using Scalar = typename DenseMatrix::template ScalarType<U>;
   public:

    struct Views
    {
      VectorView<Index> species_indices_;
      VectorView<Real> coeffs_;
      VectorView<Index> flat_ids_;

      Views() = default;

      Views(
        const Vector<Index>& species_indices,
        const Vector<Real>& coeffs,
        const Vector<Index>& flat_ids
      ) : species_indices_(species_indices.GetView()),
          coeffs_(coeffs.GetView()),
          flat_ids_(flat_ids.GetView())
      {
      }
    };
    /// @brief Name of the constraint
    std::string name_;

    /// @brief Algebraic species
    Species algebraic_species_;

    /// @brief Names of species this constraint depends on
    std::vector<std::string> species_dependencies_;

    /// @brief Species and their coefficients in the linear sum
    std::vector<StoichSpecies> terms_;

    /// @brief The constant value the linear sum should equal
    Real constant_;

    /// @brief Parameter set (unused for this class, always empty).
    std::vector<std::string> parameters_;

   protected:
    Vector<Index> species_indices_;
    Vector<Real> coeffs_;
    Vector<Index> flat_ids_;
    Views views_;

   public:
    /// @brief Default constructor
    LinearConstraint() = default;

    LinearConstraint(const LinearConstraint& other)
        : name_(other.name_),
          algebraic_species_(other.algebraic_species_),
          species_dependencies_(other.species_dependencies_),
          terms_(other.terms_),
          constant_(other.constant_),
          parameters_(other.parameters_),
          species_indices_(other.species_indices_),
          coeffs_(other.coeffs_),
          flat_ids_(other.flat_ids_),
          views_(species_indices_, coeffs_, flat_ids_)
    {
    }

    LinearConstraint& operator=(const LinearConstraint& other)
    {
      if (this != &other)
      {
        name_ = other.name_;
        algebraic_species_ = other.algebraic_species_;
        species_dependencies_ = other.species_dependencies_;
        terms_ = other.terms_;
        constant_ = other.constant_;
        parameters_ = other.parameters_;
        species_indices_ = other.species_indices_;
        coeffs_ = other.coeffs_;
        flat_ids_ = other.flat_ids_;
        views_ = Views(species_indices_, coeffs_, flat_ids_);
      }
      return *this;
    }

    LinearConstraint(LinearConstraint&& other) noexcept
      : name_(std::move(other.name_)),
        algebraic_species_(std::move(other.algebraic_species_)),
        species_dependencies_(std::move(other.species_dependencies_)),
        terms_(std::move(other.terms_)),
        constant_(std::move(other.constant_)),
        parameters_(std::move(other.parameters_)),
        species_indices_(std::move(other.species_indices_)),
        coeffs_(std::move(other.coeffs_)),
        flat_ids_(std::move(other.flat_ids_)),
        views_(species_indices_, coeffs_, flat_ids_)
    {
    }

    LinearConstraint& operator=(LinearConstraint&& other) noexcept
    {
      if (this != &other)
      {
        name_ = std::move(other.name_);
        algebraic_species_ = std::move(other.algebraic_species_);
        species_dependencies_ = std::move(other.species_dependencies_);
        terms_ = std::move(other.terms_);
        constant_ = std::move(other.constant_);
        parameters_ = std::move(other.parameters_);
        species_indices_ = std::move(other.species_indices_);
        coeffs_ = std::move(other.coeffs_);
        flat_ids_ = std::move(other.flat_ids_);
        views_ = Views(species_indices_, coeffs_, flat_ids_);
      }
      return *this;
    }

    /// @brief Construct a linear constraint
    ///        Validates that terms are non-empty
    ///        Builds species_dependencies_ from terms
    /// @param name Constraint identifier
    /// @param algebraic_species Species whose row is replaced by this algebraic constraint
    /// @param terms Vector of StoichSpecies (species, coefficient) in the linear sum
    /// @param constant The value that sum(coeff[i] * [species[i]]) should equal
    LinearConstraint(
        std::string name,
        const Species& algebraic_species,
        const std::vector<StoichSpecies>& terms,
        Real constant)
        : name_(std::move(name)),
          algebraic_species_(algebraic_species),
          terms_(terms),
          constant_(constant)
    {
      if (terms_.empty())
      {
        throw MicmException(MICM_ERROR_CATEGORY_CONSTRAINT, MICM_CONSTRAINT_ERROR_CODE_EMPTY_REACTANTS, "");
      }
      for (const auto& term : terms_)
      {
        species_dependencies_.push_back(term.species_.name_);
      }
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    /// @return Species name of the explicitly set algebraic variable
    const std::string& AlgebraicSpecies() const
    {
      return algebraic_species_.name_;
    }

    /// @brief Apply constraint parameter update (no-op for linear constraints)
    ///        Linear constraints have no temperature-dependent parameters.
    void ApplyConstraintParameter(
        const ConstraintInfo&,
        const typename DenseMatrixPolicy::template VectorType<Conditions>&,
        DenseMatrixPolicy&) const
    {
      // No-op: linear constraints don't have runtime parameters to update
    }

    void SetStateIndices(const ConstraintInfo& info, auto& jacobian_flat_ids)
    {
      std::vector<Real> coeffs_temp;
      std::vector<Index> species_indices_temp;
      for (Index i = 0; i < this->terms_.size(); ++i)
      {
        coeffs_temp.push_back(this->terms_[i].coefficient_);
        species_indices_temp.push_back(info.state_indices_[i]);
      }
      coeffs_ = coeffs_temp;
      species_indices_ = species_indices_temp;
      coeffs_.CopyToDevice();
      species_indices_.CopyToDevice();
      
      std::vector<Index> flat_ids_temp;
      flat_ids_temp.reserve(this->terms_.size());
      for (Index i = 0; i < this->terms_.size(); ++i)
      {
        flat_ids_temp.push_back(*jacobian_flat_ids++);
      }
      flat_ids_ = flat_ids_temp;
      flat_ids_.CopyToDevice();
      
      views_ = Views(species_indices_, coeffs_, flat_ids_);
    }

    /// @brief Add linear constraint residual G to forcing vector for all grid cells
    ///        Computes G = sum(coeff[i] * [species[i]]) - constant
    ///        Called directly from ConstraintSet::AddForcingTerms.
    void AddResidual(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        DenseMatrixPolicy& forcing) const
    {
      using DenseMatrix = DenseMatrixPolicy;
      using ScalarReal = typename DenseMatrix::template ScalarType<Real>;
      using ScalarIndex = typename DenseMatrix::template ScalarType<Index>;

      ScalarReal constant = this->constant_;
      ScalarIndex row_idx = info.row_index_;
      constant.CopyToDevice();
      row_idx.CopyToDevice();

      const auto& views = views_;

      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
            typename DenseMatrix::ConstViewType state_view,
            typename DenseMatrix::ConstViewType params_view,
            typename DenseMatrix::ViewType force_view)
          {
            auto linear_sum = force_view.GetRowVariable();
            state_view.ForEachRow([=](Real& sum) { sum = 0.0; }, linear_sum);

            for (Index i = 0; i < views.coeffs_.size(); ++i)
            {
              const Real coeff = views.coeffs_[i];
              const Index species_idx = views.species_indices_[i];
              state_view.ForEachRow(
                  [=](const Real& conc, Real& sum) { sum += coeff * conc; },
                  state_view.GetConstColumnView(species_idx),
                  linear_sum);
            }

            state_view.ForEachRow(
                [=](const Real& sum_val, Real& forcing_term) { forcing_term = sum_val - constant; },
                linear_sum,
                force_view.GetColumnView(row_idx));
          },
          state, state_param, forcing)(state, state_param, forcing);
    }

    /// @brief Subtract linear constraint Jacobian terms from Jacobian matrix for all grid cells
    ///        dG/d[species[i]] = coeff[i], subtracted matching SubtractJacobianTerms convention.
    ///        Called directly from ConstraintSet::SubtractJacobianTerms.
    void SubtractJacobian(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        SparseMatrixPolicy& jacobian) const
    {
      using DenseMatrix = DenseMatrixPolicy;
      using SparseMatrix = SparseMatrixPolicy;

      const auto& views = views_;

      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
            typename DenseMatrix::ConstViewType state_view,
            typename DenseMatrix::ConstViewType params_view,
            typename SparseMatrix::ViewType jacobian_values)
          {
            for (Index i = 0; i < views.coeffs_.size(); ++i)
            {
              const Real coeff = views.coeffs_[i];
              jacobian_values.ForEachBlock([=](Real& jac) { jac -= coeff; }, jacobian_values.GetBlockView(views.flat_ids_[i]));
            }
          },
          state, state_param, jacobian)(state, state_param, jacobian);
    }
  };

}  // namespace micm
