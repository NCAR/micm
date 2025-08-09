// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

// NOTE: This rate constant is a placeholder for the reaction that handles
//       the Henry's Law phase transfer. It will be refined as development moves forward.
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>

namespace micm
{
  // TODO - This is based on the implementation in cam-acom-dev
  //
  /// @brief van 't Hoff equation parameters
  /// Defines the temperature dependence of Henry's Law constants using:
  ///     H(T) = A * exp(-B * (1/T - 1/T_ref))
  /// where:
  ///     - A and B are the parameters,
  ///     - T is the current temperature,
  ///     - T_ref is the reference temperature.
  struct Van_t_HoffParams
  {
    double C = 0.0;  // mol L-1 || mol2 L-2 || mol L-1 atm-1
    double T = 0.0;  // K

    Van_t_HoffParams() = default;

    Van_t_HoffParams(double C, double T)
        : C(C),
          T(T)
    {
    }
  };

  /// @brief A rate constant for surface reactions
  class PhaseTransferRateConstant : public RateConstant
  {
   public:
    Van_t_HoffParams parameters_;

    /// @brief
    /// @param parameters The data required to build this class
    PhaseTransferRateConstant();

    /// @brief Deep copy
    std::unique_ptr<RateConstant> Clone() const override;

    std::vector<std::string> CustomParameters() const override;

    /// @brief Returns the number of custom parameters
    /// @return Number of custom parameters
    std::size_t SizeCustomParameters() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions) const override;
  };

  inline PhaseTransferRateConstant::PhaseTransferRateConstant()
  {
    // TODO
  }

  inline std::unique_ptr<RateConstant> PhaseTransferRateConstant::Clone() const
  {
    return std::make_unique<PhaseTransferRateConstant>(*this);
  }

  inline std::vector<std::string> PhaseTransferRateConstant::CustomParameters() const
  {
    return {};
  }

  inline std::size_t PhaseTransferRateConstant::SizeCustomParameters() const
  {
    return 0;
  }

  inline double PhaseTransferRateConstant::Calculate(const Conditions& conditions) const
  {
    return 0.0;
  }

  inline double PhaseTransferRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return 0.0;
  }

}  // namespace micm
