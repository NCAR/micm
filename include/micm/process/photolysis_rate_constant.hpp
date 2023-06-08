/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <string>

namespace micm
{
  class System;

  /**
   * @brief A photolysis rate constant
   */
  class PhotolysisRateConstant : public RateConstant
  {
   public:
    const std::string name_;

   public:
    /// @brief Default constructor.
    PhotolysisRateConstant();

    /// @brief
    /// @param name A name for this reaction
    PhotolysisRateConstant(const std::string& name);

    /// @brief Deep copy
    std::unique_ptr<RateConstant> clone() const override;

    /// @brief Returns the number of parameters (1) that can be set at runtime
    ///        for photolysis reactions
    ///
    ///        The single editable parameter is the unscaled rate constant for
    ///        the photolysis reaction
    /// @return Number of custom rate constant parameters
    std::size_t SizeCustomParameters() const override;

    /// @brief Calculate the rate constant
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User-defined rate constant parameters
    /// @return A rate constant based off of the conditions in the system
    double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const override;
  };

  inline PhotolysisRateConstant::PhotolysisRateConstant()
      : name_()
  {
  }

  inline PhotolysisRateConstant::PhotolysisRateConstant(const std::string& name)
      : name_(name)
  {
  }

  inline std::unique_ptr<RateConstant> PhotolysisRateConstant::clone() const
  {
    return std::unique_ptr<RateConstant>{ new PhotolysisRateConstant{ *this } };
  }

  inline double PhotolysisRateConstant::calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters)
      const
  {
    return (double)*custom_parameters;
  }

  inline std::size_t PhotolysisRateConstant::SizeCustomParameters() const
  {
    return 1;
  }
}  // namespace micm