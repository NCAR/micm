// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/phase_transfer_process.hpp>
#include <micm/process/process.hpp>
#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/yield.hpp>

#include <memory>
#include <utility>
#include <variant>
#include <vector>

namespace micm
{

  class PhaseTransferProcessBuilder
  {
   private:
    const Phase* gas_phase_;
    const Phase* condensed_phase_;
    const Phase* solvent_phase_;
    std::vector<Species> gas_species_;
    std::vector<Yield> condensed_species_;
    Species solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

   public:
    /// @brief Sets the species in the gas phase
    /// @param origin_species A vector of Species representing the origin phase species
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetGasSpecies(const Phase* phase, std::vector<Species> species)
    {
      if (!phase)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Gas Phase pointer cannot be null");

      gas_phase_ = phase;
      gas_species_ = std::move(species);
      return *this;
    }

    /// @brief Sets the species in the destination phase
    /// @param destination_species A vector of Species representing the destination phase species
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetCondensedSpecies(const Phase* phase, std::vector<Yield> condensed_species)
    {
      if (!phase)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Condensed Phase pointer cannot be null");

      condensed_phase_ = phase;
      condensed_species_ = std::move(condensed_species);
      return *this;
    }

    /// @brief Sets the solvent involved in the phase transfer process
    /// @param solvent A Species object representing the solvent
    /// @return Reference to the builder
    PhaseTransferProcessBuilder& SetSolvent(const Phase* phase, Species solvent)
    {
      if (!phase)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Solvent Phase pointer cannot be null");

      solvent_phase_ = phase;
      solvent_ = std::move(solvent);
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
      if (!coefficient_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
            "Phase Transfer Coefficient pointer cannot be null.");

      coefficient_ = std::move(coefficient);
      return *this;
    }

    /// @brief Builds the PhaseTransferProcess with the configured parameters and wraps it in a Process object
    /// @return A Process object containing the constructed PhaseTransferProcess
    Process Build()
    {
      if (!gas_phase_)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Gas Phase pointer cannot be null");
      if (!condensed_phase_)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Condensed Phase pointer cannot be null");
      if (!solvent_phase_)
        throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Solvent Phase pointer cannot be null");
      if (!coefficient_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
            "Phase Transfer Coefficient pointer cannot be null.");

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