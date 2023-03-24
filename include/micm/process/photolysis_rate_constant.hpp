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
    double rate_;
    std::string name_;

   public:
    /// @brief Default constructor.
    PhotolysisRateConstant();

    /// @brief
    /// @param rate A reaction rate, (molec cm-3)^(n-1) s-1
    PhotolysisRateConstant(double rate);

    /// @brief
    /// @param rate A reaction rate, (molec cm-3)^(n-1) s-1
    /// @param name A name for this reaction
    PhotolysisRateConstant(double rate, std::string name);

    /// @brief Calculate the rate constant
    /// @param system the system
    /// @return A rate constant based off of the conditions in the system
    double calculate(const System& system) override;
  };

  inline PhotolysisRateConstant::PhotolysisRateConstant()
      : rate_(),
        name_()
  {
  }

  inline PhotolysisRateConstant::PhotolysisRateConstant(double rate)
      : rate_(rate),
        name_()
  {
  }

  inline PhotolysisRateConstant::PhotolysisRateConstant(double rate, std::string name)
      : rate_(rate),
        name_(name)
  {
  }

  inline double PhotolysisRateConstant::calculate(const System& system)
  {
    return (double)rate_;
  }

}  // namespace micm