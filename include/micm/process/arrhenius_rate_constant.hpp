/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  /**
   * @brief An arrhenius rate constant dependent on temperature and pressure
   *
   * More information can be found here: https://open-atmos.github.io/camp/html/camp_rxn_arrhenius.html
   */
  class ArrheniusRateConstant : public RateConstant
  {
   private:
    /// @brief Pre-exponential factor, (cm−3)^(−(𝑛−1))s−1
    const double A_;
    /// @brief Unitless exponential factor
    const double B_;
    /// @brief Activation threshold, expected to be the negative activation energy divided by the boltzman constant (-E_a /
    /// k_b), K
    const double C_;
    /// @brief A factor that determines temperature dependence, (K)
    const double D_;
    /// @brief A factor that determines pressure dependence (Pa-1)
    const double E_;

   public:
    /// @brief Default constructor. All terms will be zero
    ArrheniusRateConstant();

    /// @brief An explicit constructor where each term can be set. Set B and E to zero to get the common form of the
    /// Arrhenius equation
    /// @param A Pre-exponential factor, (cm−3)^(−(𝑛−1))s−1
    /// @param B Unitless exponential factor
    /// @param C Activation threshold, expected to be the negative activation energy divided by the boltzman constant (-E_a /
    /// k_b), K
    /// @param D A factor that determines temperature dependence, (K)
    /// @param E A factor that determines pressure dependence (Pa-1)
    ArrheniusRateConstant(double A, double B, double C, double D, double E);

    /// @brief Calculate the rate constant
    /// @param system the system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const System& system) override;

    double calculate(double temperature, double pressure);
  };

  inline ArrheniusRateConstant::ArrheniusRateConstant()
      : A_(),
        B_(),
        C_(),
        D_(),
        E_()
  {
  }

  inline ArrheniusRateConstant::ArrheniusRateConstant(double A, double B, double C, double D, double E)
      : A_(A),
        B_(B),
        C_(C),
        D_(D),
        E_(E)
  {
  }

  inline double ArrheniusRateConstant::calculate(const System& system)
  {
    double temperature{ 0.1 };
    double pressure{};

    return calculate(temperature, pressure);
  }

  inline double ArrheniusRateConstant::calculate(double temperature, double pressure)
  {
    return this->A_ * std::exp(this->C_ / temperature) * pow(temperature / this->D_, this->B_) * (1.0 + this->E_ * pressure);
  }

}  // namespace micm