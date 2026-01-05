// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/transfer_coefficient/transfer_coefficient.hpp>
#include <micm/system/conditions.hpp>

#include <cmath>
#include <memory>

namespace micm
{
  /// @brief Effective Henry's Law parameters for diprotic acids
  /// For a diprotic acid like CO2: CO2(g) <-> H2CO3(aq) <-> HCO3- + H+ <-> CO32- + 2H+
  /// K_H = A * exp(C / T) = Henry's Law constant
  /// Effective Henry's Law coefficient accounts for acid dissociation:
  /// K_H_eff = K_H * (1 + K_a1/[H+] + K_a1*K_a2/[H+]^2)
  struct HenrysLawCoefficientParameters
  {
    /// @brief Pre-exponential factor for Henry's Law constant [mol L-1 atm-1]
    double A_{ 1.0 };
    /// @brief Temperature dependence factor [K]
    double C_{ 0.0 };
    /// @brief First acid dissociation constant [mol L-1]
    double K_a1_{ 0.0 };
    /// @brief Second acid dissociation constant [mol L-1]
    double K_a2_{ 0.0 };
  };

  /// @brief Effective Henry's Law coefficient for diprotic acids
  /// Accounts for temperature dependence and acid dissociation equilibria
  class HenrysLawCoefficient : public TransferCoefficient
  {
   public:
    const HenrysLawCoefficientParameters parameters_;

    /// @brief Default constructor
    HenrysLawCoefficient();

    /// @brief Explicit constructor
    /// @param parameters A set of Henry's Law coefficient parameters
    HenrysLawCoefficient(const HenrysLawCoefficientParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<TransferCoefficient> Clone() const override;

    /// @brief Calculate the effective Henry's Law coefficient
    /// @return The effective Henry's Law coefficient for the default temperature and pH
    double Calculate() const override;

    /// @brief Calculate the effective Henry's Law coefficient
    /// @param conditions The current environmental conditions of the chemical system
    /// @return The effective Henry's Law coefficient accounting for temperature and pH
    double Calculate(const Conditions& conditions) const override;

    /// @brief Calculate the effective Henry's Law coefficient
    /// @param temperature Temperature in [K]
    /// @param pH pH of the aqueous phase
    /// @return The effective Henry's Law coefficient
    double CalculateEffective(const double temperature, const double pH) const;
  };

  inline HenrysLawCoefficient::HenrysLawCoefficient()
      : parameters_()
  {
  }

  inline HenrysLawCoefficient::HenrysLawCoefficient(const HenrysLawCoefficientParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<TransferCoefficient> HenrysLawCoefficient::Clone() const
  {
    return std::make_unique<HenrysLawCoefficient>(*this);
  }

  inline double HenrysLawCoefficient::Calculate() const
  {
    // TODO (jiwon) - does this number make sense as the default?
    // Default calculation at standard temperature (298.15 K) and neutral (7.0)
    return CalculateEffective(298.15, 7.0);
  }

  inline double HenrysLawCoefficient::Calculate(const Conditions& conditions) const
  {
    if (conditions.pH.has_value())
    {
      return CalculateEffective(conditions.temperature_, conditions.pH.value());
    }
    return CalculateEffective(conditions.temperature_, 7.0);
  }

  inline double HenrysLawCoefficient::CalculateEffective(const double temperature, const double pH) const
  {
    double K_H = parameters_.A_ * std::exp(parameters_.C_ / temperature);

    // Calculate [H+] from pH
    double H_plus = std::pow(10.0, -pH);

    // Calculate effective Henry's Law coefficient for diprotic acid
    // K_H_eff = K_H * (1 + K_a1/[H+] + K_a1*K_a2/[H+]^2)
    double K_H_eff = K_H * (1.0 + parameters_.K_a1_ / H_plus + (parameters_.K_a1_ * parameters_.K_a2_) / (H_plus * H_plus));

    return K_H_eff;
  }

}  // namespace micm
