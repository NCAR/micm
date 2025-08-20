// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/species.hpp>

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
    const Phase* gas_phase_;
    const Phase* condensed_phase_;
    const Phase* solvent_phase_;
    std::vector<Species> gas_species_;
    std::vector<Yield> condensed_species_;
    Species solvent_;
    std::unique_ptr<TransferCoefficient> coefficient_;

    PhaseTransferProcess(PhaseTransferProcess&&) noexcept = default;
    PhaseTransferProcess& operator=(PhaseTransferProcess&&) noexcept = default;

    PhaseTransferProcess(
        const Phase* gas_phase,
        const Phase* condensed_phase,
        const Phase* solvent_phase,
        std::vector<Species> gas_species,
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
        if (!other.gas_phase_)
          throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Gas Phase pointer cannot be null");
        if (!other.condensed_phase_)
          throw std::system_error(
              make_error_code(MicmProcessErrc::PhaseIsNotSet), "Condensed Phase pointer cannot be null");
        if (!other.solvent_phase_)
          throw std::system_error(make_error_code(MicmProcessErrc::PhaseIsNotSet), "Solvent Phase pointer cannot be null");
        if (!other.coefficient_)
          throw std::system_error(
              make_error_code(MicmProcessErrc::TransferCoefficientIsNotSet),
              "Phase Transfer Coefficient pointer cannot be null.");

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
    }
  };

}  // namespace micm
