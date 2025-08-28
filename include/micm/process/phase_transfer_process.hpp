// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_error.hpp>
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

  /// @brief Represents a phase-transfer process where species moves from one phase to another
  class PhaseTransferProcess
  {
   public:
    Phase gas_phase_;
    Phase condensed_phase_;
    Phase solvent_phase_;
    Species gas_species_;
    std::vector<Yield> condensed_species_;
    Species solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

    PhaseTransferProcess(PhaseTransferProcess&&) noexcept = default;
    PhaseTransferProcess& operator=(PhaseTransferProcess&&) noexcept = default;

    PhaseTransferProcess(
        const Phase& gas_phase,
        const Phase& condensed_phase,
        const Phase& solvent_phase,
        Species gas_species,
        std::vector<Yield> condensed_species,
        Species solvent,
        std::unique_ptr<TransferCoefficient> coefficient)
        : gas_phase_(gas_phase),
          condensed_phase_(condensed_phase),
          solvent_phase_(solvent_phase),
          gas_species_(std::move(gas_species)),
          condensed_species_(std::move(condensed_species)),
          solvent_(std::move(solvent)),
          coefficient_(std::move(coefficient))
    {
      Validate();
    }

    PhaseTransferProcess(const PhaseTransferProcess& other)
        : gas_phase_(other.gas_phase_),
          condensed_phase_(other.condensed_phase_),
          solvent_phase_(other.solvent_phase_),
          gas_species_(other.gas_species_),
          condensed_species_(other.condensed_species_),
          solvent_(other.solvent_),
          coefficient_(other.coefficient_ ? other.coefficient_->Clone() : nullptr)
    {
      Validate();
    }

    PhaseTransferProcess& operator=(const PhaseTransferProcess& other)
    {
      if (this != &other)
      {
        if (!other.coefficient_)
          throw std::system_error(
              make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
              "Cannot copy from a PhaseTransferProcess with null coefficient");

        gas_phase_ = other.gas_phase_;
        condensed_phase_ = other.condensed_phase_;
        solvent_phase_ = other.solvent_phase_;
        gas_species_ = other.gas_species_;
        condensed_species_ = other.condensed_species_;
        solvent_ = other.solvent_;
        coefficient_ = other.coefficient_->Clone();
      }

      return *this;
    }

   private:
    void Validate() const
    {
      if (!coefficient_)
        throw std::system_error(
            make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
            "Phase Transfer Coefficient pointer cannot be null");
    }
  };

}  // namespace micm
