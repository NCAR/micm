/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  struct TroeRateConstantParameters
  {
    /// @brief // TODO:
    double k0_A_;
    /// @brief // TODO:
    double k0_B_;
    /// @brief // TODO:
    double k0_C_;
    /// @brief // TODO:
    double kinf_A_;
    /// @brief // TODO:
    double kinf_B_;
    /// @brief // TODO:
    double kinf_C_;
    /// @brief // TODO:
    double Fc_;
    /// @brief // TODO:
    double N_;
  };

  /**
   * @brief A Troe rate constant
   *
   */
  class TroeRateConstant : public RateConstant
  {
   public:
    const TroeRateConstantParameters parameters_;

   public:
    /// @brief Default constructor is not allowed
    TroeRateConstant() = delete;

    /// @brief An explicit constructor
    /// @param parameters A set of troe rate constants
    TroeRateConstant(TroeRateConstantParameters parameters);

    /// @brief Calculate the rate constant
    /// @param system the system
    /// @return A rate constant based off of the conditions in the systßem
    double calculate(const System& system) override;

    /// @brief Calculate the rate constant
    /// @param temperature Temperature in [K]
    /// @param air_number_density Number density in [# cm-3]
    /// @return
    double calculate(double temperature, double air_number_density);
  };

  inline TroeRateConstant::TroeRateConstant(TroeRateConstantParameters parameters)
      : parameters_(parameters)
  {
  }

  inline double TroeRateConstant::calculate(const System& system)
  {
    double temperature{}, air_number_density{};

    return calculate(temperature, air_number_density);
  }

  inline double TroeRateConstant::calculate(double temperature, double air_number_density)
  {
    double k0 = parameters_.k0_A_ * std::exp(parameters_.k0_C_ / temperature) * pow(temperature / 300.0, parameters_.k0_B_);
    double kinf =
        parameters_.kinf_A_ * std::exp(parameters_.kinf_C_ / temperature) * pow(temperature / 300.0, parameters_.kinf_B_);

    return k0 * air_number_density / (1.0 + k0 * air_number_density / kinf) *
           pow(parameters_.Fc_, 1.0 / (1.0 + 1.0 / parameters_.N_ * pow(log10(k0 * air_number_density / kinf), 2)));
  }

}  // namespace micm