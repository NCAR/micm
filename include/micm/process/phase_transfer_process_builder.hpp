// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/phase_transfer_process.hpp>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/process/process.hpp>
#include <micm/system/species.hpp>

#include <memory>
#include <utility>
#include <vector>
#include <variant>

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
    PhaseTransferProcessBuilder& SetOriginSpecies(const std::vector<SpeciesInPhase>& origin_species)
    {
      origin_species_ = origin_species;
      return *this;
    }

    PhaseTransferProcessBuilder& SetDestinationSpecies(const std::vector<SpeciesInPhase>& destination_species)
    {
      destination_species_ = destination_species;
      return *this;
    }

    PhaseTransferProcessBuilder& SetSolvent(const SpeciesInPhase& solvent)
    {
      solvent_ = solvent;
      return *this;
    }

    PhaseTransferProcessBuilder& SetTransferCoefficient(const TransferCoefficient& coefficient)
    {
      coefficient_ = coefficient.Clone();
      return *this;
    }

    PhaseTransferProcessBuilder& SetTransferCoefficient(std::unique_ptr<TransferCoefficient> coefficient)
    {
      coefficient_ = std::move(coefficient);
      return *this;
    }

    /// @brief Build the phase transfer process and wrap in Process
    Process Build()
    {
      PhaseTransferProcess process(std::move(origin_species_),
                                   std::move(destination_species_),
                                   std::move(solvent_),
                                   std::move(coefficient_));
      return Process(std::move(process));
    }

  };

}  // namespace micm