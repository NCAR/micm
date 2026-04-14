// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/process.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class ChemicalReactionBuilder
  {
   public:
    /// @brief Sets the list of reactant species involved in the chemical reaction.
    /// @param reactants A list of Species objects representing the reactants
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetReactants(const std::vector<Species>& reactants)
    {
      reactants_ = reactants;
      return *this;
    }

    /// @brief Sets the list of product species and their yields for the chemical reaction.
    /// @param products A list of StoichSpecies objects representing the products
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetProducts(const std::vector<StoichSpecies>& products)
    {
      products_ = products;
      return *this;
    }

    /// @brief Sets the rate constant from any supported parameter struct.
    ///        Accepts any type that is a member of RateConstantVariant.
    /// @param rate_constant Parameter struct (e.g. ArrheniusRateConstantParameters)
    /// @return Reference to the builder
    template<class T>
    ChemicalReactionBuilder& SetRateConstant(T&& rate_constant)
    {
      rate_constant_ = RateConstantVariant(std::forward<T>(rate_constant));
      return *this;
    }

    /// @brief Sets the phase in which the chemical reaction occurs (e.g., gas, aqueous)
    /// @param phase Phase object representing the reaction phase
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetPhase(const Phase& phase)
    {
      phase_ = phase;
      return *this;
    }

    /// @brief Transfers ownership of all internally stored data into a ChemicalReaction,
    ///        then wraps it into a Process using std::variant
    /// @return A Process containing the constructed ChemicalReaction
    /// @throws MicmException if the rate constant has not been set
    Process Build()
    {
      if (!rate_constant_.has_value())
        throw MicmException(
            MicmSeverity::Error,
            MICM_ERROR_CATEGORY_PROCESS,
            MICM_PROCESS_ERROR_CODE_RATE_CONSTANT_IS_NOT_SET,
            "Rate constant has not been set.");

      ChemicalReaction reaction(
          std::move(reactants_), std::move(products_), std::move(rate_constant_.value()), phase_);
      return Process(std::move(reaction));
    }

   private:
    std::vector<Species> reactants_;
    std::vector<StoichSpecies> products_;
    std::optional<RateConstantVariant> rate_constant_;
    Phase phase_;
  };

}  // namespace micm
