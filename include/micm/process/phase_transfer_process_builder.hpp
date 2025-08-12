// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/phase_transfer_process.hpp>
#include <micm/process/process.hpp>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/species.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class PhaseTransferProcessBuilder
  {
   private:
    std::vector<SpeciesInPhase> origin_species_;
    std::vector<SpeciesInPhase> destination_species_;
    SpeciesInPhase solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

   public:
    /// @brief Sets the species in the origin phase
    /// @param origin_species A vector of SpeciesInPhase representing the origin phase species
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetOriginSpecies(const std::vector<SpeciesInPhase>& origin_species)
    {
      origin_species_ = origin_species;
      return *this;
    }

    /// @brief Sets the species in the destination phase
    /// @param destination_species A vector of SpeciesInPhase representing the destination phase species
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetDestinationSpecies(const std::vector<SpeciesInPhase>& destination_species)
    {
      destination_species_ = destination_species;
      return *this;
    }

    /// @brief Sets the solvent involved in the phase transfer process
    /// @param solvent A SpeciesInPhase object representing the solvent
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetSolvent(const SpeciesInPhase& solvent)
    {
      solvent_ = solvent;
      return *this;
    }

    /// @brief Sets the transfer coefficient by cloning the provided coefficient object
    /// @param coefficient A reference to a TransferCoefficient to be cloned
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetTransferCoefficient(const TransferCoefficient& coefficient)
    {
      coefficient_ = coefficient.Clone();
      return *this;
    }

    /// @brief Sets the transfer coefficient by taking ownership of a unique_ptr
    ///        This method avoids cloning and assumes ownership of the given pointer
    /// @param coefficient A unique_ptr to a TransferCoefficient object
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetTransferCoefficient(std::unique_ptr<TransferCoefficient> coefficient)
    {
      coefficient_ = std::move(coefficient);
      return *this;
    }

    /// @brief Builds the PhaseTransferProcess with the configured parameters and wraps it in a Process object
    /// @return A Process object containing the constructed PhaseTransferProcess
    Process Build()
    {
      PhaseTransferProcess process(
          std::move(origin_species_), std::move(destination_species_), std::move(solvent_), std::move(coefficient_));
      return Process(std::move(process));
    }
  };

}  // namespace micm