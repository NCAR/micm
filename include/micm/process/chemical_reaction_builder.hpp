// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/process.hpp>
#include <micm/system/species.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/yield.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class ChemicalReactionBuilder
  {
   private:
    Phase phase_;
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;

   public:
    /// @brief Sets the phase name for the chemical reaction (e.g., "gas", "aqueous")
    /// @param phase A string representing the phase in which the reaction occurs
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetPhase(const Phase& phase)
    {
      phase_ = phase;
      return *this;
    }

    /// @brief Sets the list of reactant species involved in the chemical reaction
    /// @param reactants A vector of Species objects representing the reactants
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetReactants(const std::vector<Species>& reactants)
    {
      reactants_ = reactants;
      return *this;
    }

    /// @brief Sets the list of product species and their yields for the chemical reaction
    /// @param products A vector of Yield objects representing the products
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetProducts(const std::vector<Yield>& products)
    {
      products_ = products;
      return *this;
    }

    /// @brief Sets the rate constant by cloning the provided RateConstant object
    ///        This method performs a deep copy of the given rate constant using its Clone() method.
    ///        Useful when the original rate constant must remain unchanged.
    /// @param rate_constant A reference to a RateConstant object to be cloned
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetRateConstant(const RateConstant& rate_constant)
    {
      rate_constant_ = rate_constant.Clone();
      return *this;
    }

    ///@brief Sets the rate constant by transferring ownership of a RateConstant object
    ///       This method avoids cloning and takes ownership of the provided unique_ptr.
    ///@param rate_constant A unique_ptr to a RateConstant
    ///@return Reference to the builder
    ChemicalReactionBuilder& SetRateConstant(std::unique_ptr<RateConstant> rate_constant)
    {
      rate_constant_ = std::move(rate_constant);
      return *this;
    }

    /// @brief Transfers ownership of all internally stored data into a ChemicalReaction,
    ///        then wraps it into a Process using std::variant
    /// @return A Process containing the constructed ChemicalReaction
    Process Build()
    {
      ChemicalReaction reaction(
          std::move(phase_), std::move(reactants_), std::move(products_), std::move(rate_constant_));
      return Process(std::move(reaction));
    }
  };

}  // namespace micm