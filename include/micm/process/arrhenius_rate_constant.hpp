/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>

namespace micm
{

  /**
   * @brief An arrhenius rate constant
   *
   * @tparam T The type of each factor of the arrhenius rate constant
   */
  template<typename T>
  class ArrheniusRateConstant : public RateConstant
  {
   private:
    /// @brief //TODO
    const T A_;
    /// @brief //TODO
    const T B_;
    /// @brief //TODO
    const T C_;
    /// @brief //TODO
    const T D_;
    /// @brief //TODO
    const T E_;

   public:
    /// @brief Default constructor
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