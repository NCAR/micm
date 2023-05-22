/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <utility>
#include <vector>

namespace micm
{

  using Yield = std::pair<micm::Species, double>;

  Yield yields(micm::Species species, double yield)
  {
    return Yield(species, yield);
  };

  struct Process
  {
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    RateConstant rate_constant_;
    Phase* phase_;
  };
}  // namespace micm
