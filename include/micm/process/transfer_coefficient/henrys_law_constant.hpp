// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/constants.hpp>

#include <cmath>
#include <memory>

namespace micm
{

  /// @brief Henry's Law constant parameters
  struct HenrysLawConstantParameters
  {
    /// @brief Henry's Law constant for the species at the reference temperature [mol L-1 atm-1]
    double H_ref_{ 1.3e-3 };
    /// @brief Enthalpy of dissolution for the species in solution [J/mol]
    double enthalpy_{ -12000.0 };
    /// @brief The standard reference temperature [K]
    double temperature_ref_ {298.15};
  };

  /// @brief Henry's Law constant that accounts for temperature dependence
  class HenrysLawConstant : public TransferCoefficient
  {
   public:
    const HenrysLawConstantParameters parameters_;

    /// @brief Default constructor
    HenrysLawConstant()
      : parameters_()
    {
    }
    
    /// @brief Explicit constructor
    /// @param parameters A set of Henry's Law constant parameters
    HenrysLawConstant(const HenrysLawConstantParameters& parameters)
      : parameters_(parameters)
    {
    }

    /// @brief Deep copy
    std::unique_ptr<TransferCoefficient> Clone() const override
    {
      return std::make_unique<HenrysLawConstant>(*this);
    }

    /// @brief Calculate the Henry's Law constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return The Henry's Law constant accounting for temperature
    double Calculate(const Conditions& conditions) const override
    {
      return parameters_.H_ref_ * std::exp((-parameters_.enthalpy_ / constants::GAS_CONSTANT)
              * (1.0 / conditions.temperature_ - 1.0 / parameters_.temperature_ref_));
    }

  };

}  // namespace micm


