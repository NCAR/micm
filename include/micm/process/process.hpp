/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <vector>
#include <micm/system/species.hpp>
#include <micm/system/phase.hpp>
#include <micm/process/rate_constant.hpp>

namespace micm
{
  struct Process {
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    RateConstant rate_constant_;
    Phase* phase_;
  };
}  // namespace micm
