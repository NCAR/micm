/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

namespace micm
{

  class System;

  /**
   * @brief A base class for any type of rate constant
   *
   */
  class RateConstant
  {
   public:
    /// @brief Virtual destructor
    virtual ~RateConstant(){};
    /// @brief Calculate the rate constant for a set of conditions
    virtual double calculate(const System& system) { return 0; };

   private:
  };

}  // namespace micm