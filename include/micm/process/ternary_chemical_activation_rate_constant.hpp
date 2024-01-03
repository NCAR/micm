// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  struct TernaryChemicalActivationRateConstantParameters
  {
    /// @brief low-pressure pre-exponential factor
    double k0_A_ = 1.0;
    /// @brief low-pressure temperature-scaling parameter
    double k0_B_ = 0.0;
    /// @brief low-pressure exponential factor
    double k0_C_ = 0.0;
    /// @brief high-pressure pre-exponential factor
    double kinf_A_ = 1.0;
    /// @brief high-pressure temperature-scaling parameter
    double kinf_B_ = 0.0;
    /// @brief high-pressure exponential factor
    double kinf_C_ = 0.0;
    /// @brief TernaryChemicalActivation F_c parameter
    double Fc_ = 0.6;
    /// @brief TernaryChemicalActivation N parameter
    double N_ = 1.0;
  };

  /// @brief A TernaryChemicalActivation rate constant
  class TernaryChemicalActivationRateConstant : public RateConstant
  {
   public:
    const TernaryChemicalActivationRateConstantParameters parameters_;

    /// @brief Default constructor
    TernaryChemicalActivationRateConstant();

    /// @brief An explicit constructor
    /// @param parameters A set of troe rate constants
    TernaryChemicalActivationRateConstant(const TernaryChemicalActivationRateConstantParameters& parameters);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions) const override;

    /// @brief Calculate the rate constant
    /// @param temperature Temperature in [K]
    /// @param air_number_density Number density in [mol m-3]
    /// @return
    double calculate(const double& temperature, const double& air_number_density) const;
  };

  inline TernaryChemicalActivationRateConstant::TernaryChemicalActivationRateConstant()
      : parameters_()
  {
  }

  inline TernaryChemicalActivationRateConstant::TernaryChemicalActivationRateConstant(
      const TernaryChemicalActivationRateConstantParameters& parameters)
      : parameters_(parameters)
  {
  }

  inline std::unique_ptr<RateConstant> TernaryChemicalActivationRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new TernaryChemicalActivationRateConstant{ *this } };
  }

  inline double TernaryChemicalActivationRateConstant::calculate(const Conditions& conditions) const
  {
    double val = calculate(conditions.temperature_, conditions.air_density_);
    return val;
  }

  inline double TernaryChemicalActivationRateConstant::calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    double val = calculate(conditions.temperature_, conditions.air_density_);
    return val;
  }

  inline double TernaryChemicalActivationRateConstant::calculate(const double& temperature, const double& air_number_density)
      const
  {
    double k0 =
        parameters_.k0_A_ * std::exp(parameters_.k0_C_ / temperature) * std::pow(temperature / 300.0, parameters_.k0_B_);
    double kinf = parameters_.kinf_A_ * std::exp(parameters_.kinf_C_ / temperature) *
                  std::pow(temperature / 300.0, parameters_.kinf_B_);

    return k0 / (1.0 + k0 * air_number_density / kinf) *
           std::pow(
               parameters_.Fc_, parameters_.N_ / (parameters_.N_ + std::pow(std::log10(k0 * air_number_density / kinf), 2)));
  }

}  // namespace micm
