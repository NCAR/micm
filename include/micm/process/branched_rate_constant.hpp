// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <micm/process/rate_constant.hpp>
#include <micm/util/constants.hpp>

namespace micm
{

  struct BranchedRateConstantParameters
  {
    enum class Branch
    {
      Alkoxy,
      Nitrate
    };
    /// @brief reaction branch
    Branch branch_;
    /// @brief pre-exponential factor
    double X_;
    /// @brief exponential factor
    double Y_;
    /// @brief branching factor
    double a0_;
    /// @brief number of heavy atoms in the RO2 reacting species (excluding the peroxy moiety)
    int n_;
  };

  /// @brief A Branched rate constant
  class BranchedRateConstant : public RateConstant
  {
   public:
    const BranchedRateConstantParameters parameters_;
    const double k0_;
    const double z_;

    /// @brief Default constructor
    BranchedRateConstant();

    /// @brief An explicit constructor
    /// @param parameters A set of branched rate constant parameters
    BranchedRateConstant(const BranchedRateConstantParameters& parameters);

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

    /// @brief Calculate A(T,[M],n)
    /// @param temperature Temperature in [K]
    /// @param air_number_density Number density of air in [mol m-3]
    double A(const double& temperature, const double& air_number_density) const;
  };

  inline BranchedRateConstant::BranchedRateConstant()
      : parameters_(),
        k0_(),
        z_()
  {
  }

  inline BranchedRateConstant::BranchedRateConstant(const BranchedRateConstantParameters& parameters)
      : parameters_(parameters),
        k0_(2.0e-22 * AVOGADRO_CONSTANT * 1.0e-6 * std::exp(parameters_.n_)),
        z_(A(293.0, 2.45e19 / AVOGADRO_CONSTANT * 1.0e6) * (1.0 - parameters_.a0_) / parameters_.a0_)
  {
  }

  inline std::unique_ptr<RateConstant> BranchedRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new BranchedRateConstant{ *this } };
  }

  inline double BranchedRateConstant::calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    return calculate(conditions.temperature_, conditions.air_density_);
  }

  inline double BranchedRateConstant::calculate(const Conditions& conditions) const
  {
    return calculate(conditions.temperature_, conditions.air_density_);
  }

  inline double BranchedRateConstant::calculate(const double& temperature, const double& air_number_density) const
  {
    double pre = parameters_.X_ * std::exp(-parameters_.Y_ / temperature);
    double Atmn = A(temperature, air_number_density);
    if (parameters_.branch_ == BranchedRateConstantParameters::Branch::Alkoxy)
      return pre * (z_ / (z_ + Atmn));
    return pre * (Atmn / (Atmn + z_));
  }

  inline double BranchedRateConstant::A(const double& temperature, const double& air_number_density) const
  {
    double a = k0_ * air_number_density;
    double b = 0.43 * std::pow(temperature / 298.0, -8);
    return a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2)));
  }
}  // namespace micm