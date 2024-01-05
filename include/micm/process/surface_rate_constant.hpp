// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include <cmath>
#include <micm/process/rate_constant.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/property_keys.hpp>
#include <string>

namespace micm
{
  struct SurfaceRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Gas-phase species reacting on surface
    Species species_;
    /// @brief Reaction probability (0-1) [unitless]
    double reaction_probability_{ 1.0 };
  };

  /// @brief A rate constant for surface reactions
  class SurfaceRateConstant : public RateConstant
  {
   public:
    SurfaceRateConstantParameters parameters_;
    double diffusion_coefficient_;   // [m2 s-1]
    double mean_free_speed_factor_;  // 8 * gas_constant / ( pi * molecular_weight )  [K-1]

    /// @brief
    /// @param parameters The data required to build this class
    SurfaceRateConstant(const SurfaceRateConstantParameters& parameters);

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
    double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions) const override;
  };

  inline SurfaceRateConstant::SurfaceRateConstant(const SurfaceRateConstantParameters& parameters)
      : parameters_(parameters),
        diffusion_coefficient_(parameters.species_.properties_.at(GAS_DIFFUSION_COEFFICIENT)),
        mean_free_speed_factor_(8.0 * GAS_CONSTANT / (M_PI * parameters.species_.properties_.at(MOLECULAR_WEIGHT)))
  {
  }

  inline std::unique_ptr<RateConstant> SurfaceRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new SurfaceRateConstant{ *this } };
  }

  inline double SurfaceRateConstant::calculate(const Conditions& conditions) const
  {
    throw std::runtime_error(
        "Surface rate constants must be supplied with a radius and number density using the alternative calculate function");
  }

  inline double SurfaceRateConstant::calculate(
      const Conditions& conditions,
      std::vector<double>::const_iterator custom_parameters) const
  {
    const double mean_free_speed = std::sqrt(mean_free_speed_factor_ * conditions.temperature_);
    const double radius = *(custom_parameters++);
    const double number = *(custom_parameters);
    double val = (double)4.0 * number * M_PI * radius * radius /
                 (radius / diffusion_coefficient_ + 4.0 / (mean_free_speed * parameters_.reaction_probability_));
    return val;
  }

  inline std::vector<std::string> SurfaceRateConstant::CustomParameters() const
  {
    return std::vector<std::string>{ parameters_.label_ + ".effective radius [m]",
                                     parameters_.label_ + ".particle number concentration [# m-3]" };
  }

  inline std::size_t SurfaceRateConstant::SizeCustomParameters() const
  {
    return 2;
  }
}  // namespace micm
