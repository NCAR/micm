/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  class IntraPhaseProcess
  {
   private:
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    RateConstant rate_constant_;

   public:
  };

}  // namespace micm