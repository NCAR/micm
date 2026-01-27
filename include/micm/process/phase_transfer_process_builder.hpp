// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/phase_transfer_process.hpp>
#include <micm/process/process.hpp>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class PhaseTransferProcessBuilder
  {
   private:
    Phase gas_phase_;
    Phase condensed_phase_;
    Phase solvent_phase_;
    Species gas_species_;
    Species condensed_species_;
    Species solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

   public:
    /// @brief Sets the species in the gas phase
    /// @param phase Phase object representing the gas phase
    /// @param species Species object in the gas phase
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetGasSpecies(const Phase& phase, const Species& species)
    {
      gas_phase_ = phase;
      gas_species_ = species;
      return *this;
    }

    /// @brief Sets the species in the condensed phase
    /// @param phase Phase object representing the condensed phase
    /// @param condensed_species Species object in the condensed phase
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetCondensedSpecies(const Phase& phase, const Species& species)
    {
      condensed_phase_ = phase;
      condensed_species_ = species;
      return *this;
    }

    /// @brief Sets the solvent involved in the phase transfer process
    /// @param phase Phase object representing the solvent phase
    /// @param solvent A Species object representing the solvent
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetSolvent(const Phase& phase, const Species& solvent)
    {
      solvent_phase_ = phase;
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

    /// @brief Builds the PhaseTransferProcess with the configured parameters and wraps it in a Process object
    /// @return A Process object containing the constructed PhaseTransferProcess
    /// @throws std::system_error if the provided coefficient pointer is null
    Process Build()
    {
      if (!coefficient_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
            "Phase Transfer Coefficient pointer cannot be null");

      PhaseTransferProcess process(
          gas_phase_,
          condensed_phase_,
          solvent_phase_,
          std::move(gas_species_),
          std::move(condensed_species_),
          std::move(solvent_),
          std::move(coefficient_));
      return Process(std::move(process));
    }
  };

}  // namespace micm