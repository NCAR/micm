// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <memory>
#include <vector>

namespace micm
{
  /// @brief Represents a species in a chemical reaction, defined by its stoichiometric coefficient
  struct StoichSpecies
  {
    Species species_;
    double coefficient_{ 1.0 };
  };

}  // namespace micm