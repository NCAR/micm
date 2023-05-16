/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/condition.hpp>

namespace micm
{

  /**
   * @brief A rate constant from an external model
   *
   */
  class ExternalRateConstant : public RateConstant
  {
   private:
    /// @brief The rate
    const double rate_;
    /// @brief The condition this rate applies to
    const Condition condition_;

   public:
    /// @brief Default constructor
    ExternalRateConstant();

    /// @brief Calculate a reaction rate
    /// @param system A defined system
    /// @return A reaction rate
    double calculate(const System& system) override;
  };

  inline ExternalRateConstant::ExternalRateConstant()
      : rate_(),
        condition_("", "")
  {
  }

  inline double ExternalRateConstant::calculate(const System& system)
  {
    return 0.0;
  }

}  // namespace micm