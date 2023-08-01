// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <math.h>
#include <micm/process/rate_constant.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/property_keys.hpp>
#include <string>

namespace micm
{
  class System;

  /// @brief A rate constant for surface reactions
  class SurfaceRateConstant : public RateConstant
  {
    std::string name_;
    double diffusion_coefficient_;   // [m2 s-1]
    double mean_free_speed_factor_;  // 8 * gas_constant / ( pi * molecular_weight )  [K-1]
    double reaction_probability_;    // reaction probability (0-1)

   public:
    /// @brief Default constructor.
    SurfaceRateConstant();

    /// @brief
    /// @param name A name for this reaction
    SurfaceRateConstant(const std::string& name, const Species& species, const double reaction_probability);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

    /// @brief Returns labels for the surface rate constant parameters:
    ///        aerosol effective radius [m]
    ///        aerosol number concentration [# m-3]
    /// @return Rate constant labels
    std::vector<std::string> CustomParameters() const override;

    /// @brief Returns the number of custom parameters
    /// @return Number of custom parameters
    std::size_t SizeCustomParameters() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters)
        const override;
  };

  inline SurfaceRateConstant::SurfaceRateConstant()
      : name_(),
        diffusion_coefficient_(),
        mean_free_speed_factor_(),
        reaction_probability_()
  {
  }

  inline SurfaceRateConstant::SurfaceRateConstant(
      const std::string& name,
      const Species& species,
      const double reaction_probability = 1.0)
      : name_(name),
        diffusion_coefficient_(species.properties_.at(GAS_DIFFUSION_COEFFICIENT)),
        mean_free_speed_factor_(8.0 * GAS_CONSTANT / (M_PI * species.properties_.at(MOLECULAR_WEIGHT))),
        reaction_probability_(reaction_probability)
  {
  }

  inline std::unique_ptr<RateConstant> SurfaceRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new SurfaceRateConstant{ *this } };
  }

  inline double SurfaceRateConstant::calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    const double mean_free_speed = std::sqrt(mean_free_speed_factor_ * conditions.temperature_);
    const double radius = *(custom_parameters++);
    const double number = *(custom_parameters);
    return (double)4.0 * number * M_PI * radius * radius /
           (radius / diffusion_coefficient_ + 4.0 / (mean_free_speed * reaction_probability_));
  }

  inline std::vector<std::string> SurfaceRateConstant::CustomParameters() const
  {
    return std::vector<std::string>{ "SURF." + name_ + ".effective radius [m]",
                                     "SURF." + name_ + ".particle number concentration [# m-3]" };
  }

  inline std::size_t SurfaceRateConstant::SizeCustomParameters() const
  {
    return 2;
  }
}  // namespace micm