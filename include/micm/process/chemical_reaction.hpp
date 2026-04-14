// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/taylor_series_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>

#include <variant>
#include <vector>

namespace micm
{

  /// @brief Value-typed union of all supported rate constant parameter types.
  ///        Stored by value in ChemicalReaction; consumed at store-build time by
  ///        ReactionRateStore::BuildFrom.  Never sent to GPU.
  using RateConstantVariant = std::variant<
      ArrheniusRateConstantParameters,
      TroeRateConstantParameters,
      TernaryChemicalActivationRateConstantParameters,
      BranchedRateConstantParameters,
      TunnelingRateConstantParameters,
      TaylorSeriesRateConstantParameters,
      ReversibleRateConstantParameters,
      UserDefinedRateConstantParameters,
      SurfaceRateConstantParameters,
      LambdaRateConstantParameters>;

  /// @brief Represents a chemical reaction with reactants, products, rate constant and phase
  class ChemicalReaction
  {
   public:
    std::vector<Species> reactants_;
    std::vector<StoichSpecies> products_;
    RateConstantVariant rate_constant_;
    Phase phase_;

    ChemicalReaction() = default;
    ChemicalReaction(const ChemicalReaction&) = default;
    ChemicalReaction(ChemicalReaction&&) noexcept = default;
    ChemicalReaction& operator=(const ChemicalReaction&) = default;
    ChemicalReaction& operator=(ChemicalReaction&&) noexcept = default;

    ChemicalReaction(
        std::vector<Species> reactants,
        std::vector<StoichSpecies> products,
        RateConstantVariant rate_constant,
        const Phase& phase)
        : reactants_(std::move(reactants)),
          products_(std::move(products)),
          rate_constant_(std::move(rate_constant)),
          phase_(phase)
    {
    }
  };

}  // namespace micm
