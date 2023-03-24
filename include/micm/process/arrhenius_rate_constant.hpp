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
    ArrheniusRateConstantParameters parameters_;

   public:
    /// @brief Default constructor. All terms will be zero
    ArrheniusRateConstant();

    /// @brief An explicit constructor where each term can be set. Set B and E to zero to get the common form of the
    /// Arrhenius equation
    /// @param parameters A set of arrhenius rate constants
    ArrheniusRateConstant(ArrheniusRateConstantParameters parameters);

    /// @brief Copy constructor
    /// @param other
    ArrheniusRateConstant(const ArrheniusRateConstant& other);

    /// @brief Move constructor
    /// @param other
    ArrheniusRateConstant(ArrheniusRateConstant&& other);

    /// @brief Copy assignment operator
    ArrheniusRateConstant& operator=(ArrheniusRateConstant& other);

    /// @brief Move assignment operator
    ArrheniusRateConstant& operator=(ArrheniusRateConstant&& other);

    /// @brief Calculate the rate constant
    /// @param system the system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const System& system) override;

    double calculate(double temperature, double pressure);
  };

  inline ArrheniusRateConstant::ArrheniusRateConstant()
      : parameters_()
  {
  }

  inline ArrheniusRateConstant::ArrheniusRateConstant(ArrheniusRateConstantParameters parameters)
      : parameters_(parameters)
  {
  }

  inline ArrheniusRateConstant::ArrheniusRateConstant(const ArrheniusRateConstant& other)
      : parameters_(other.parameters_)
  {
  }

  inline ArrheniusRateConstant::ArrheniusRateConstant(ArrheniusRateConstant&& other)
      : parameters_(std::move(other.parameters_))
  {
  }

  inline ArrheniusRateConstant& ArrheniusRateConstant::operator=(ArrheniusRateConstant& other)
  {
    if (this != &other)
    {
      parameters_ = other.parameters_;
    }
    return *this;
  }

  inline ArrheniusRateConstant& ArrheniusRateConstant::operator=(ArrheniusRateConstant&& other)
  {
    if (this != &other)
    {
      parameters_ = std::move(other.parameters_);
    }
    return *this;
  }

  inline double ArrheniusRateConstant::calculate(const System& system)
  {
    double temperature{ 0.1 };
    double pressure{};

    return calculate(temperature, pressure);
  }

  inline double ArrheniusRateConstant::calculate(double temperature, double pressure)
  {
    return parameters_.A_ * std::exp(parameters_.C_ / temperature) * pow(temperature / parameters_.D_, parameters_.B_) *
           (1.0 + parameters_.E_ * pressure);
  }

}  // namespace micm