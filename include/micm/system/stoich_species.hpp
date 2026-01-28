// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/species.hpp>

namespace micm
{
  /// @brief Represents a species in a chemical reaction, defined by its stoichiometric coefficient
  struct StoichSpecies
  {
    Species species_;
    double coefficient_{ 1.0 };
  };

  [[deprecated("micm::Yield has been renamed to micm::StoichSpecies; please use StoichSpecies instead")]]
  using Yield = StoichSpecies;

}  // namespace micm