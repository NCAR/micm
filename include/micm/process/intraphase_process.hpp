/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  /**
   * @brief An intraphase process
   *
   */
  template<class Rate>
  class IntraPhaseProcess
  {
   private:
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    Rate rate_constant_;

   public:
  };

}  // namespace micm