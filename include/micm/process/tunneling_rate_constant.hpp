// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant.hpp>

#include <cmath>

namespace micm
{

  struct TunnelingRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    double A_ = 1.0;
    /// @brief Linear temperature-dependent parameter [K]
    double B_ = 0.0;
    /// @brief Cubed temperature-dependent parameter [K^3]
    double C_ = 0.0;
  };

  /// @brief Rate constant for tunneling reactions
  class TunnelingRateConstant : public RateConstant
  {
   public:
    const TunnelingRateConstantParameters parameters_;

    /// @brief Default constructor
    TunnelingRateConstant();

    /// @brief An explicit constructor
    /// @param parameters A set of troe rate constants
    TunnelingRateConstant(const TunnelingRateConstantParameters& parameters);

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
    /// @return the calculated rate constant
    double Calculate(const double& temperature) const;
  };

  inline TunnelingRateConstant::TunnelingRateConstant()
      : parameters_()
  {
  }

  inline TunnelingRateConstant::TunnelingRateConstant(const TunnelingRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> TunnelingRateConstant::Clone() const
  {
    return std::unique_ptr<RateConstant>{ new TunnelingRateConstant{ *this } };
  }

  inline double TunnelingRateConstant::Calculate(const Conditions& conditions) const
  {
    return Calculate(conditions.temperature_);
  }

  inline double TunnelingRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return Calculate(conditions.temperature_);
  }

  inline double TunnelingRateConstant::Calculate(const double& temperature) const
  {
    return parameters_.A_ * std::exp(-parameters_.B_ / temperature + parameters_.C_ / std::pow(temperature, 3));
  }

}  // namespace micm