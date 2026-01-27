// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/chemical_reaction.hpp>
#include <micm/process/process.hpp>
#include <micm/process/process_error.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/utils.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class ChemicalReactionBuilder
  {
   public:
    /// @brief Enables aerosol scoping for reactant and product species
    ///        This function must be called before setting reactants or products
    ///        in order for scoping to be applied.
    ///        Cannot be used together with SetPhase. They are mutually exclusive.
    /// @param scope Aerosol scope prefix to apply to species names
    /// @param phase Phase object representing the reaction phase
    /// @return Reference to the builder
    /// @throws std::system_error if SetPhase, SetReactants, or SetProducts has already been called
    ChemicalReactionBuilder& SetAerosolScope(const std::string& scope, const Phase& phase)
    {
      if (has_phase_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::InvalidConfiguration),
            "SetPhase and SetAerosolScope are mutually exclusive and should not be used together.");
      if (has_reactants_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::InvalidConfiguration), "SetAerosolScope must be called before SetReactants.");
      if (has_products_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::InvalidConfiguration), "SetAerosolScope must be called before SetProducts.");

      scope_ = scope;
      phase_ = phase;
      has_scope_ = true;

      return *this;
    }

    /// @brief Sets the list of reactant species involved in the chemical reaction.
    ///        When scoping is enabled, each reactant name is prefixed with the preset scope.
    /// @param reactants A list of Species objects representing the reactants
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetReactants(const std::vector<Species>& reactants)
    {
      if (has_scope_)
      {
        reactants_.reserve(reactants.size());
        for (const auto& species : reactants)
        {
          reactants_.push_back(species);
          Scope(reactants_.back(), phase_);
        }
      }
      else
      {
        reactants_ = reactants;
      }

      has_reactants_ = true;

      return *this;
    }

    /// @brief Sets the list of product species and their yields for the chemical reaction.
    ///        When scoping is enabled, each product name is prefixed with the preset scope.
    /// @param products A list of StoichSpecies objects representing the products
    /// @return Reference to the builder
    ChemicalReactionBuilder& SetProducts(const std::vector<StoichSpecies>& products)
    {
      if (has_scope_)
      {
        products_.reserve(products.size());
        for (const auto& [species, coefficient] : products)
        {
          products_.emplace_back(species, coefficient);
          Scope(products_.back().species_, phase_);
        }
      }
      else
      {
        products_ = products;
      }

      has_products_ = true;

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

    /// @brief Sets the phase in which the chemical reaction occurs (e.g., gas, aqueous)
    ///        Cannot be used together with SetAerosolScope. They are mutually exclusive.
    /// @param phase Phase object representing the reaction phase
    /// @return Reference to the builder
    /// @throws std::system_error if SetAerosolScope has already been called
    ChemicalReactionBuilder& SetPhase(const Phase& phase)
    {
      if (has_scope_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::InvalidConfiguration),
            "SetPhase and SetAerosolScope are mutually exclusive and should not be used together.");

      phase_ = phase;
      has_phase_ = true;
      return *this;
    }

    /// @brief Transfers ownership of all internally stored data into a ChemicalReaction,
    ///        then wraps it into a Process using std::variant
    /// @return A Process containing the constructed ChemicalReaction
    /// @throws std::system_error if the provided rate constant pointer is null
    Process Build()
    {
      if (!rate_constant_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::RateConstantIsNotSet), "Rate Constant pointer cannot be null.");

      ChemicalReaction reaction(std::move(reactants_), std::move(products_), std::move(rate_constant_), phase_);
      return Process(std::move(reaction));
    }

   private:
    std::vector<Species> reactants_;
    std::vector<StoichSpecies> products_;
    std::unique_ptr<RateConstant> rate_constant_;
    Phase phase_;
    std::string scope_;

    bool has_scope_ = false;
    bool has_phase_ = false;
    bool has_reactants_ = false;
    bool has_products_ = false;

    /// @brief Applies an aerosol phase-specific scope to a species by prefixing its name
    /// @param species Species object whose name will be modified
    /// @param phase Phase whose name is used in the scope prefix
    void Scope(Species& species, const Phase& phase)
    {
      species.name_ = JoinStrings({ scope_, phase.name_, species.name_ });
    }
  };

}  // namespace micm