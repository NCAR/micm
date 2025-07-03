// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant.hpp>

#include <cmath>

namespace micm
{
  struct TaylorSeriesRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    double A_{ 1 };
    /// @brief Unitless exponential factor
    double B_{ 0 };
    /// @brief Activation threshold, expected to be the negative activation energy divided by the boltzman constant
    ///        [-E_a / k_b), K]
    double C_{ 0 };
    /// @brief A factor that determines temperature dependence [K]
    double D_{ 300 };
    /// @brief A factor that determines pressure dependence [Pa-1]
    double E_{ 0 };
    /// @brief Taylor coefficients for the series expansion
    std::vector<double> coefficients_{ 1.0 };
  };

  /// @brief An taylor series rate constant dependent on temperature and pressure
  class TaylorSeriesRateConstant : public RateConstant
  {
   public:
    const TaylorSeriesRateConstantParameters parameters_;

    /// @brief Default constructor
    TaylorSeriesRateConstant();

    /// @brief An explicit constructor where each term can be set. Set B and E to zero to get the common form of the
    /// TaylorSeries equation
    /// @param parameters A set of arrhenius rate constants
    TaylorSeriesRateConstant(const TaylorSeriesRateConstantParameters& parameters);

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

    double Calculate(const double& temperature, const double& pressure) const;
  };

  inline TaylorSeriesRateConstant::TaylorSeriesRateConstant()
      : parameters_()
  {
  }

  inline TaylorSeriesRateConstant::TaylorSeriesRateConstant(const TaylorSeriesRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> TaylorSeriesRateConstant::Clone() const
  {
    return std::unique_ptr<RateConstant>{ new TaylorSeriesRateConstant{ *this } };
  }

  inline double TaylorSeriesRateConstant::Calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return Calculate(conditions.temperature_, conditions.pressure_);
  }

  inline double TaylorSeriesRateConstant::Calculate(const Conditions& conditions) const
  {
    return Calculate(conditions.temperature_, conditions.pressure_);
  }

  inline double TaylorSeriesRateConstant::Calculate(const double& temperature, const double& pressure) const
  {
    double result = 0.0;
    for (size_t i = 0; i < parameters_.coefficients_.size(); ++i)
    {
      result += parameters_.coefficients_[i] * std::pow(temperature, i);
    }
    return result * parameters_.A_ * std::exp(parameters_.C_ / temperature) * std::pow(temperature / parameters_.D_, parameters_.B_) *
           (1.0 + parameters_.E_ * pressure);
  }

}  // namespace micm