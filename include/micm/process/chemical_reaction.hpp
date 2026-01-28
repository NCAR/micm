// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_error.hpp>
#include <micm/process/rate_constant/rate_constant.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/error.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  /// @brief Represents a chemical reaction with reactants, products, rate constant and phase
  class ChemicalReaction
  {
   public:
    std::vector<Species> reactants_;
    std::vector<StoichSpecies> products_;
    std::unique_ptr<RateConstant> rate_constant_;
    Phase phase_;

    ChemicalReaction(ChemicalReaction&&) noexcept = default;
    ChemicalReaction& operator=(ChemicalReaction&&) noexcept = default;

    ChemicalReaction(
        std::vector<Species> reactants,
        std::vector<StoichSpecies> products,
        std::unique_ptr<RateConstant> rate_constant,
        const Phase& phase)
        : reactants_(std::move(reactants)),
          products_(std::move(products)),
          rate_constant_(std::move(rate_constant)),
          phase_(phase)
    {
      Validate();
    }

    ChemicalReaction(const ChemicalReaction& other)
        : reactants_(other.reactants_),
          products_(other.products_),
          rate_constant_(other.rate_constant_ ? other.rate_constant_->Clone() : nullptr),
          phase_(other.phase_)
    {
      Validate();
    }

    ChemicalReaction& operator=(const ChemicalReaction& other)
    {
      if (this != &other)
      {
        if (!other.rate_constant_)
          throw std::system_error(
              make_error_code(MicmProcessErrc::RateConstantIsNotSet),
              "Cannot copy from a ChemicalReaction with null rate constant");

        reactants_ = other.reactants_;
        products_ = other.products_;
        rate_constant_ = other.rate_constant_->Clone();
        phase_ = other.phase_;
      }

      return *this;
    }

    /// @brief Calculates the rate constants for each process for the current state
    /// @param processes The set of processes for which rate constants are to be calculated
    /// @param state The current solver state that will be modified with the updated rate constants
    /// This function is overloaded based on whether DenseMatrixPolicy is vectorizable.
    template<
        class DenseMatrixPolicy,
        class SparseMatrixPolicy,
        class LuDecompositionPolicy,
        class LMatrixPolicy,
        class UMatrixPolicy>
      requires(!VectorizableDense<DenseMatrixPolicy>)
    static void CalculateRateConstants(
        const std::vector<ChemicalReaction>& processes,
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state);

    template<
        class DenseMatrixPolicy,
        class SparseMatrixPolicy,
        class LuDecompositionPolicy,
        class LMatrixPolicy,
        class UMatrixPolicy>
      requires(VectorizableDense<DenseMatrixPolicy>)
    static void CalculateRateConstants(
        const std::vector<ChemicalReaction>& processes,
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state);

   private:
    void Validate() const
    {
      if (!rate_constant_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::RateConstantIsNotSet), "Rate Constant pointer cannot be null");
    }
  };

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
    requires(!VectorizableDense<DenseMatrixPolicy>)
  void ChemicalReaction::CalculateRateConstants(
      const std::vector<ChemicalReaction>& processes,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state)
  {
    for (std::size_t i{}; i < state.custom_rate_parameters_.NumRows(); ++i)
    {
      const std::vector<double> custom_parameters = state.custom_rate_parameters_[i];
      std::vector<double>::const_iterator custom_parameters_iter = custom_parameters.begin();
      std::size_t i_rate_constant = 0;
      for (auto& process : processes)
      {
        double fixed_reactants = 1.0;
        for (auto& reactant : process.reactants_)
          if (reactant.IsParameterized())
            fixed_reactants *= reactant.parameterize_(state.conditions_[i]);
        state.rate_constants_[i][(i_rate_constant++)] =
            process.rate_constant_->Calculate(state.conditions_[i], custom_parameters_iter) * fixed_reactants;
        custom_parameters_iter += process.rate_constant_->SizeCustomParameters();
      }
    }
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
    requires(VectorizableDense<DenseMatrixPolicy>)
  void ChemicalReaction::CalculateRateConstants(
      const std::vector<ChemicalReaction>& processes,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state)
  {
    const auto& v_custom_parameters = state.custom_rate_parameters_.AsVector();
    auto& v_rate_constants = state.rate_constants_.AsVector();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    // loop over all rows
    for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
    {
      std::size_t offset_rc = i_group * state.rate_constants_.GroupSize();
      std::size_t offset_params = i_group * state.custom_rate_parameters_.GroupSize();
      std::size_t rate_const_size = std::min(L, state.rate_constants_.NumRows() - (i_group * L));
      for (auto& process : processes)
      {
        std::vector<double> params(process.rate_constant_->SizeCustomParameters());
        for (std::size_t i_cell{}; i_cell < rate_const_size; ++i_cell)
        {
          for (std::size_t i_param = 0; i_param < params.size(); ++i_param)
          {
            params[i_param] = v_custom_parameters[offset_params + i_param * L + i_cell];
          }
          std::vector<double>::const_iterator custom_parameters_iter = params.begin();
          double fixed_reactants = 1.0;
          for (auto& reactant : process.reactants_)
            if (reactant.IsParameterized())
              fixed_reactants *= reactant.parameterize_(state.conditions_[i_group * L + i_cell]);
          v_rate_constants[offset_rc + i_cell] =
              process.rate_constant_->Calculate(state.conditions_[i_group * L + i_cell], custom_parameters_iter) *
              fixed_reactants;
        }
        offset_params += params.size() * L;
        offset_rc += L;
      }
    }
  }

}  // namespace micm