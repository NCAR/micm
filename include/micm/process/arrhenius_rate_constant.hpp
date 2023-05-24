/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  class System;

  struct ArrheniusRateConstantParameters
  {
    /// @brief Pre-exponential factor, (cm−3)^(−(𝑛−1))s−1
    double A_{ 1 };
    /// @brief Unitless exponential factor
    double B_{ 0 };
    /// @brief Activation threshold, expected to be the negative activation energy divided by the boltzman constant (-E_a /
    /// k_b), K
    double C_{ 0 };
    /// @brief A factor that determines temperature dependence, (K)
    double D_{ 300 };
    /// @brief A factor that determines pressure dependence (Pa-1)
    double E_{ 0 };
  };

  /**
   * @brief An arrhenius rate constant dependent on temperature and pressure
   *
   * More information can be found here: https://open-atmos.github.io/camp/html/camp_rxn_arrhenius.html
   */
  class ArrheniusRateConstant : public RateConstant
  {
   public:
    const ArrheniusRateConstantParameters parameters_;

   public:
    /// @brief Default constructor. All terms will be zero
    ArrheniusRateConstant();

    /// @brief An explicit constructor where each term can be set. Set B and E to zero to get the common form of the
    /// Arrhenius equation
    /// @param parameters A set of arrhenius rate constants
    ArrheniusRateConstant(ArrheniusRateConstantParameters parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

    /// @brief Calculate the rate constant
    /// @param state The current state of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double calculate(const State& state, std::vector<double>::const_iterator custom_parameters) const override;

    double calculate(const double& temperature, const double& pressure) const;
  };

  inline ArrheniusRateConstant::ArrheniusRateConstant()
      : parameters_()
  {
  }

  inline ArrheniusRateConstant::ArrheniusRateConstant(ArrheniusRateConstantParameters parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> ArrheniusRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new ArrheniusRateConstant{ *this } };
  }

  inline double ArrheniusRateConstant::calculate(const State& state, std::vector<double>::const_iterator custom_parameters)
      const
  {
    return calculate(state.temperature_, state.pressure_);
  }

  inline double ArrheniusRateConstant::calculate(const double& temperature, const double& pressure) const
  {
    return parameters_.A_ * std::exp(parameters_.C_ / temperature) * pow(temperature / parameters_.D_, parameters_.B_) *
           (1.0 + parameters_.E_ * pressure);
  }

}  // namespace micm