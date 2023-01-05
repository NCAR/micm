/* Copyright (C) 2022 National Center for Atmospheric Research,
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
   * @tparam T The type of the rate constant
   */
  template<typename T>
  class ExternalRateConstant : public RateConstant
  {
   private:
    /// @brief The rate
    const T rate_;
    /// @brief The condition this rate applies to
    const Condition condition_;

   public:
    /// @brief Default constructor
    ExternalRateConstant();
  };

  template<typename T>
  inline ExternalRateConstant<T>::ExternalRateConstant()
      : rate_(),
        condition_()
  {
  }

}  // namespace micm