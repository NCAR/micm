// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/species.hpp>

#include <memory>
#include <utility>
#include <vector>
#include <variant>


namespace micm
{

  struct SpeciesInPhase
  {
    std::string phase_name_;
    Species species_;

    SpeciesInPhase() = default;

    SpeciesInPhase(const std::string& phase_name, const Species& species)
        : phase_name_(phase_name),
          species_(species)
    {
    }

  };

  /// @brief Represents a phase-transfer process where reactants move from one phase to another,
  ///        producing products in the destination phase
  class PhaseTransferProcess
  {
  public:
    std::vector<SpeciesInPhase> origin_species_;
    std::vector<SpeciesInPhase> destination_species_;
    SpeciesInPhase solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

    PhaseTransferProcess(PhaseTransferProcess&&) noexcept = default;
    PhaseTransferProcess& operator=(PhaseTransferProcess&&) noexcept = default;

    PhaseTransferProcess(std::vector<SpeciesInPhase> origin_species,
                         std::vector<SpeciesInPhase> destination_species,
                         SpeciesInPhase solvent,
                         std::unique_ptr<TransferCoefficient> coefficient)
      : origin_species_(std::move(origin_species)),
        destination_species_(std::move(destination_species)),
        solvent_(std::move(solvent)),
        coefficient_(std::move(coefficient))
    {}

    PhaseTransferProcess(const PhaseTransferProcess& other)
      : origin_species_(other.origin_species_),
        destination_species_(other.destination_species_),
        solvent_(other.solvent_),
        coefficient_(other.coefficient_ ? other.coefficient_->Clone() : nullptr)
    {}

    PhaseTransferProcess& operator=(const PhaseTransferProcess& other)
    {
      if (this != &other)
      {
        origin_species_ = other.origin_species_;
        destination_species_ = other.destination_species_;
        solvent_ = other.solvent_;
        coefficient_ = other.coefficient_ ? other.coefficient_->Clone() : nullptr;
      }

      return *this;
    }

  };

}  // namespace micm
