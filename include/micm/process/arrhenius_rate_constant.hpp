/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>

namespace micm
{

  template<typename T>
  class ArrheniusRateConstant : public RateConstant
  {
   private:
    const T A_;
    const T B_;
    const T C_;
    const T D_;
    const T E_;

   public:
    ArrheniusRateConstant();
  };

  template<typename T>
  inline ArrheniusRateConstant<T>::ArrheniusRateConstant()
      : A_(),
        B_(),
        C_(),
        D_(),
        E_()
  {
  }

}  // namespace micm