// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint_info.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <system_error>
#include <vector>

namespace micm
{

  /// @brief Constraint for chemical equilibrium with temperature-dependent K_eq using Van't Hoff equation
  ///        For a reversible reaction: aA + bB <-> cC + dD
  ///        The equilibrium constraint is: G = K_eq(T) * [A]^a * [B]^b - [C]^c * [D]^d = 0
  ///        where K_eq(T) = K_HLC_ref * exp((delta_H / R) * (1/T - 1/T_ref))
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class EquilibriumConstraint
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
    /// @brief Define parameters for Van't Hoff equation
    struct VantHoffParam
    {
      Real K_HLC_ref_;                    // Henry’s Law constant at the reference temperature (typically 298.15K)
      Real delta_H_;                      // Enthalpy of dissolution (J/mol)
      Real R_ = constants::GAS_CONSTANT;  // (J/mol·K)
      Real T_ref_ = 298.15;
    };

    struct Views
    {
      VectorView<Real> reactant_stoich_;
      VectorView<Index> reactant_state_idx_;
      VectorView<Real> product_stoich_;
      VectorView<Index> product_state_idx_;
      VectorView<Index> flat_ids_;

      Views() = default;

      Views(
          const Vector<Real>& reactant_stoich,
          const Vector<Index>& reactant_state_idx,
          const Vector<Real>& product_stoich,
          const Vector<Index>& product_state_idx,
          const Vector<Index>& flat_ids)
          : reactant_stoich_(reactant_stoich.GetView()),
            reactant_state_idx_(reactant_state_idx.GetView()),
            product_stoich_(product_stoich.GetView()),
            product_state_idx_(product_state_idx.GetView()),
            flat_ids_(flat_ids.GetView())
      {
      }
    };

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

   private:
    /// @brief Van't Hoff equation parameter used to calculate Henry’s Law Constant
    VantHoffParam vant_hoff_param_;

    /// @brief Indices into the reactants_ vector for each species dependency
    std::vector<Index> reactant_dependency_indices_;

    /// @brief Indices into the products_ vector for each species dependency
    std::vector<Index> product_dependency_indices_;

   protected:
    Vector<Real> reactant_stoich_;
    Vector<Index> reactant_state_idx_;
    Vector<Real> product_stoich_;
    Vector<Index> product_state_idx_;
    Vector<Index> flat_ids_;
    Views views_;

   public:
    /// @brief Default constructor
    EquilibriumConstraint() = default;

    EquilibriumConstraint(const EquilibriumConstraint& other)
        : name_(other.name_),
          algebraic_species_(other.algebraic_species_),
          species_dependencies_(other.species_dependencies_),
          reactants_(other.reactants_),
          products_(other.products_),
          parameters_(other.parameters_),
          vant_hoff_param_(other.vant_hoff_param_),
          reactant_dependency_indices_(other.reactant_dependency_indices_),
          product_dependency_indices_(other.product_dependency_indices_),
          reactant_stoich_(other.reactant_stoich_),
          reactant_state_idx_(other.reactant_state_idx_),
          product_stoich_(other.product_stoich_),
          product_state_idx_(other.product_state_idx_),
          flat_ids_(other.flat_ids_),
          views_(reactant_stoich_, reactant_state_idx_, product_stoich_, product_state_idx_, flat_ids_)
    {
    }

    EquilibriumConstraint& operator=(const EquilibriumConstraint& other)
    {
      if (this != &other)
      {
        name_ = other.name_;
        algebraic_species_ = other.algebraic_species_;
        species_dependencies_ = other.species_dependencies_;
        reactants_ = other.reactants_;
        products_ = other.products_;
        parameters_ = other.parameters_;
        vant_hoff_param_ = other.vant_hoff_param_;
        reactant_dependency_indices_ = other.reactant_dependency_indices_;
        product_dependency_indices_ = other.product_dependency_indices_;
        reactant_stoich_ = other.reactant_stoich_;
        reactant_state_idx_ = other.reactant_state_idx_;
        product_stoich_ = other.product_stoich_;
        product_state_idx_ = other.product_state_idx_;
        flat_ids_ = other.flat_ids_;
        views_ = Views(reactant_stoich_, reactant_state_idx_, product_stoich_, product_state_idx_, flat_ids_);
      }
      return *this;
    }

    EquilibriumConstraint(EquilibriumConstraint&& other) noexcept
        : name_(std::move(other.name_)),
          algebraic_species_(other.algebraic_species_),  // Species is not move constructible
          species_dependencies_(std::move(other.species_dependencies_)),
          reactants_(std::move(other.reactants_)),
          products_(std::move(other.products_)),
          parameters_(std::move(other.parameters_)),
          vant_hoff_param_(std::move(other.vant_hoff_param_)),
          reactant_dependency_indices_(std::move(other.reactant_dependency_indices_)),
          product_dependency_indices_(std::move(other.product_dependency_indices_)),
          reactant_stoich_(std::move(other.reactant_stoich_)),
          reactant_state_idx_(std::move(other.reactant_state_idx_)),
          product_stoich_(std::move(other.product_stoich_)),
          product_state_idx_(std::move(other.product_state_idx_)),
          flat_ids_(std::move(other.flat_ids_)),
          views_(reactant_stoich_, reactant_state_idx_, product_stoich_, product_state_idx_, flat_ids_)
    {
    }

    EquilibriumConstraint& operator=(EquilibriumConstraint&& other) noexcept
    {
      if (this != &other)
      {
        name_ = std::move(other.name_);
        algebraic_species_ = other.algebraic_species_;  // Speces is not move assignable
        species_dependencies_ = std::move(other.species_dependencies_);
        reactants_ = std::move(other.reactants_);
        products_ = std::move(other.products_);
        parameters_ = std::move(other.parameters_);
        vant_hoff_param_ = std::move(other.vant_hoff_param_);
        reactant_dependency_indices_ = std::move(other.reactant_dependency_indices_);
        product_dependency_indices_ = std::move(other.product_dependency_indices_);
        reactant_stoich_ = std::move(other.reactant_stoich_);
        reactant_state_idx_ = std::move(other.reactant_state_idx_);
        product_stoich_ = std::move(other.product_stoich_);
        product_state_idx_ = std::move(other.product_state_idx_);
        flat_ids_ = std::move(other.flat_ids_);
        views_ = Views(reactant_stoich_, reactant_state_idx_, product_stoich_, product_state_idx_, flat_ids_);
      }
      return *this;
    }

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
          vant_hoff_param_(std::move(vant_hoff_param))
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
    }

    /// @brief Returns the species whose row should be replaced by this algebraic constraint
    /// @return Species name of the explicitly set algebraic variable
    const std::string& AlgebraicSpecies() const
    {
      return algebraic_species_.name_;
    }

    /// @brief Apply temperature-dependent K_eq parameter update for each grid cell
    ///        Computes K_eq(T) using the Van't Hoff equation and writes to state_param[K_eq_idx].
    ///        Called directly from ConstraintSet::UpdateStateParameters before each solve.
    /// @param info Constraint information including state parameter indices
    /// @param conditions Per-grid-cell atmospheric conditions (temperature, pressure, etc.)
    /// @param state_param State parameter matrix to update
    void ApplyConstraintParameter(
        const ConstraintInfo& info,
        const typename DenseMatrixPolicy::template VectorType<Conditions>& conditions,
        DenseMatrixPolicy& state_param) const
    {
      Index K_eq_idx = info.state_param_indices_[0];
      const VantHoffParam& param = vant_hoff_param_;
      
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename Vector<Conditions>::ConstViewType& conditions_view,
              const typename DenseMatrixPolicy::ViewType& state_param_view) {
            state_param_view.ForEachRow(
                [=](const Conditions& cond, Real& K_eq)
                { K_eq = param.K_HLC_ref_ * std::exp((param.delta_H_ / param.R_) * (1.0 / cond.temperature_ - 1.0 / param.T_ref_)); },
                conditions_view,
                state_param_view.GetColumnView(K_eq_idx));
          },
          conditions,
          state_param)(conditions, state_param);
    }

    void SetStateIndices(const ConstraintInfo& info, auto& jacobian_flat_ids)
    {
      std::vector<Index> flat_ids_temp;
      flat_ids_temp.reserve(reactants_.size() + products_.size());
      for (Index i = 0; i < reactants_.size(); ++i)
      {
        flat_ids_temp.push_back(*jacobian_flat_ids++);
      }
      for (Index i = 0; i < products_.size(); ++i)
      {
        flat_ids_temp.push_back(*jacobian_flat_ids++);
      }
      flat_ids_ = flat_ids_temp;
      flat_ids_.CopyToDevice();

      std::vector<Real> reactant_stoich_temp;
      std::vector<Index> reactant_state_idx_temp;
      for (Index i = 0; i < this->reactants_.size(); ++i)
      {
        reactant_stoich_temp.push_back(this->reactants_[i].coefficient_);
        reactant_state_idx_temp.push_back(info.state_indices_[reactant_dependency_indices_[i]]);
      }
      reactant_stoich_ = reactant_stoich_temp;
      reactant_state_idx_ = reactant_state_idx_temp;
      reactant_stoich_.CopyToDevice();
      reactant_state_idx_.CopyToDevice();

      std::vector<Real> product_stoich_temp;
      std::vector<Index> product_state_idx_temp;
      for (Index i = 0; i < this->products_.size(); ++i)
      {
        product_stoich_temp.push_back(this->products_[i].coefficient_);
        product_state_idx_temp.push_back(info.state_indices_[product_dependency_indices_[i]]);
      }
      product_stoich_ = product_stoich_temp;
      product_state_idx_ = product_state_idx_temp;
      product_stoich_.CopyToDevice();
      product_state_idx_.CopyToDevice();

      views_ = Views(reactant_stoich_, reactant_state_idx_, product_stoich_, product_state_idx_, flat_ids_);
    }

    /// @brief Create function object to compute equilibrium constraint residual for all grid cells
    ///        Computes G = K_eq(T) * prod([reactants]^stoich) - prod([products]^stoich) for the algebraic constraint
    ///        Called during solver build (SetConstraintFunctions) to pre-compile residual computation
    /// @param info Constraint information including row index, species indices, and parameter indices
    /// @brief Add equilibrium constraint residual G to forcing vector for all grid cells
    ///        Computes G = K_eq(T) * prod([reactants]^stoich) - prod([products]^stoich)
    ///        Called directly from ConstraintSet::AddForcingTerms.
    /// @param info Constraint information including row index and parameter indices
    /// @param state Current species concentrations
    /// @param state_param Current state parameters (contains K_eq column)
    /// @param forcing Forcing terms — constraint row is overwritten with residual G
    void AddResidual(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        DenseMatrixPolicy& forcing) const
    {
      Index row_idx = info.row_index_;
      Index K_eq_idx = info.state_param_indices_[0];
      
      const auto& views = views_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& state_param_view,
              const typename DenseMatrixPolicy::ViewType& force_view) {
            auto reactant_product = force_view.GetRowVariable();
            auto product_product = force_view.GetRowVariable();

            state_view.ForEachRow(
                [=](const Real& K_eq, Real& rp, Real& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                state_param_view.GetConstColumnView(K_eq_idx),
                reactant_product,
                product_product);

            for (Index i = 0; i < views.reactant_stoich_.size(); ++i)
            {
              const Real stoich = views.reactant_stoich_[i];
              const Index species_idx = views.reactant_state_idx_[i];
              state_view.ForEachRow(
                  [=](const Real& conc, Real& product) { product *= std::pow((conc > Real(0) ? conc : Real(0)), stoich); },
                  state_view.GetConstColumnView(species_idx),
                  reactant_product);
            }

            for (Index i = 0; i < views.product_stoich_.size(); ++i)
            {
              const Real stoich = views.product_stoich_[i];
              const Index species_idx = views.product_state_idx_[i];
              state_view.ForEachRow(
                  [=](const Real& conc, Real& product) { product *= std::pow((conc > Real(0) ? conc : Real(0)), stoich); },
                  state_view.GetConstColumnView(species_idx),
                  product_product);
            }

            state_view.ForEachRow(
                [=](const Real& rp, const Real& pp, Real& forcing_term) { forcing_term = rp - pp; },
                reactant_product,
                product_product,
                force_view.GetColumnView(row_idx));
          },
          state,
          state_param,
          forcing)(state, state_param, forcing);
    }

    /// @brief Subtract Jacobian partial derivatives dG/d[species] from Jacobian matrix for all grid cells
    ///        Called directly from ConstraintSet::SubtractJacobianTerms.
    /// @param info Constraint information including row index and parameter indices
    /// @param state Current species concentrations
    /// @param state_param Current state parameters (contains K_eq column)
    /// @param jacobian Sparse Jacobian matrix to update
    void SubtractJacobian(
        const ConstraintInfo& info,
        const DenseMatrixPolicy& state,
        const DenseMatrixPolicy& state_param,
        SparseMatrixPolicy& jacobian) const
    {
      Index K_eq_idx = info.state_param_indices_[0];

      const auto& views = views_;

      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrix::ConstViewType& state_view,
              const typename DenseMatrix::ConstViewType& state_param_view,
              const typename SparseMatrix::ViewType& jacobian_values) {
            auto reactant_product = jacobian_values.GetBlockVariable();
            auto product_product = jacobian_values.GetBlockVariable();
            auto partial_derivative = jacobian_values.GetBlockVariable();

            jacobian_values.ForEachBlock(
                [=](const Real& K_eq, Real& rp, Real& pp)
                {
                  rp = K_eq;
                  pp = 1.0;
                },
                state_param_view.GetConstColumnView(K_eq_idx),
                reactant_product,
                product_product);

            for (Index i = 0; i < views.reactant_stoich_.size(); ++i)
            {
              const Real stoich = views.reactant_stoich_[i];
              const Index species_idx = views.reactant_state_idx_[i];
              jacobian_values.ForEachBlock(
                  [=](const Real& conc, Real& product) { product *= std::pow((conc > Real(0) ? conc : Real(0)), stoich); },
                  state_view.GetConstColumnView(species_idx),
                  reactant_product);
            }

            for (Index i = 0; i < views.product_stoich_.size(); ++i)
            {
              const Real stoich = views.product_stoich_[i];
              const Index species_idx = views.product_state_idx_[i];
              jacobian_values.ForEachBlock(
                  [=](const Real& conc, Real& product) { product *= std::pow((conc > Real(0) ? conc : Real(0)), stoich); },
                  state_view.GetConstColumnView(species_idx),
                  product_product);
            }

            for (Index i = 0; i < views.reactant_stoich_.size(); ++i)
            {
              const Real stoich_i = views.reactant_stoich_[i];
              const Index species_idx_i = views.reactant_state_idx_[i];

              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock(
                  [=](const Real& K_eq, Real& prod) { prod = K_eq; },
                  state_param_view.GetConstColumnView(K_eq_idx),
                  partial_product);

              for (Index j = 0; j < views.reactant_stoich_.size(); ++j)
              {
                if (j != i)
                {
                  const Real stoich_j = views.reactant_stoich_[j];
                  const Index species_idx_j = views.reactant_state_idx_[j];
                  jacobian_values.ForEachBlock(
                      [=](const Real& conc, Real& prod) { prod *= std::pow((conc > Real(0) ? conc : Real(0)), stoich_j); },
                      state_view.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              jacobian_values.ForEachBlock(
                  [=](const Real& conc, const Real& prod, Real& partial)
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
                  state_view.GetConstColumnView(species_idx_i),
                  partial_product,
                  partial_derivative);

              jacobian_values.ForEachBlock(
                  [=](const Real& partial, Real& jac) { jac -= partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(views.flat_ids_[i]));
            }

            for (Index i = 0; i < views.product_stoich_.size(); ++i)
            {
              const Real stoich_i = views.product_stoich_[i];
              const Index species_idx_i = views.product_state_idx_[i];

              auto partial_product = jacobian_values.GetBlockVariable();
              jacobian_values.ForEachBlock([=](Real& prod) { prod = 1.0; }, partial_product);

              for (Index j = 0; j < views.product_stoich_.size(); ++j)
              {
                if (j != i)
                {
                  const Real stoich_j = views.product_stoich_[j];
                  const Index species_idx_j = views.product_state_idx_[j];
                  jacobian_values.ForEachBlock(
                      [=](const Real& conc, Real& prod) { prod *= std::pow((conc > Real(0) ? conc : Real(0)), stoich_j); },
                      state_view.GetConstColumnView(species_idx_j),
                      partial_product);
                }
              }

              jacobian_values.ForEachBlock(
                  [=](const Real& conc, const Real& prod, Real& partial)
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
                  state_view.GetConstColumnView(species_idx_i),
                  partial_product,
                  partial_derivative);

              jacobian_values.ForEachBlock(
                  [=](const Real& partial, Real& jac) { jac += partial; },
                  partial_derivative,
                  jacobian_values.GetBlockView(views.flat_ids_[views.reactant_stoich_.size() + i]));
            }
          },
          state,
          state_param,
          jacobian)(state, state_param, jacobian);
    }
  };

}  // namespace micm
