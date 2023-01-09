/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/condition.hpp>

namespace micm
{

  template<typename T>
  class ExternalRateConstant : public RateConstant
  {
   private:
    const T rate_;
    const Condition condition_;

  public:
    ExternalRateConstant();
  };

  template<typename T>
  inline ExternalRateConstant<T>::ExternalRateConstant()
    : rate_(), condition_()
  {
  }

}  // namespace micm