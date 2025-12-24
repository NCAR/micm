// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant/rate_constant.hpp>

#include <cmath>

namespace micm
{
  struct ReversibleRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    double A_{ 1 };
    /// @brief Activation threshold [K]
    double C_{ 0 };
    /// @brief Reverse rate constant [s−1]
    double k_r_{ 0 };
  };

  /// @brief A reversible rate constant with temperature dependence
  class ReversibleRateConstant : public RateConstant
  {
   public:
    const ReversibleRateConstantParameters parameters_;

    /// @brief Default constructor
    ReversibleRateConstant();

    /// @brief An explicit constructor
    /// @param parameters A set of reversible rate constant parameters
    ReversibleRateConstant(const ReversibleRateConstantParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> Clone() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double Calculate(const Conditions& conditions) const override;

    /// @brief Calculate the rate constant
    /// @param temperature Temperature in [K]
    /// @return A rate constant based on temperature
    double Calculate(const double& temperature) const;
  };

  inline ReversibleRateConstant::ReversibleRateConstant()
      : parameters_()
  {
  }

  inline ReversibleRateConstant::ReversibleRateConstant(const ReversibleRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> ReversibleRateConstant::Clone() const
  {
    return std::make_unique<ReversibleRateConstant>(*this);
  }

  inline double ReversibleRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return Calculate(conditions.temperature_);
  }

  inline double ReversibleRateConstant::Calculate(const Conditions& conditions) const
  {
    return Calculate(conditions.temperature_);
  }

  // Condensed phase reversible reaction
  // K_eq = A * exp(C / T) = Equilibrium constant
  // k_r = reverse rate constant
  // (k_f = K_eq * k_r = forward rate constant)
  inline double ReversibleRateConstant::Calculate(const double& temperature) const
  {
    double K_eq = parameters_.A_ * std::exp(parameters_.C_ / temperature);
    return K_eq * parameters_.k_r_;
  }

}  // namespace micm
