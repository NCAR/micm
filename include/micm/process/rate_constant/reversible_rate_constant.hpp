// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

namespace micm
{
  struct ReversibleRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    Real A_{ 1 };
    /// @brief Activation threshold [K]
    Real C_{ 0 };
    /// @brief Reverse rate constant [s−1], indicating how fast the species
    ///        leaves the condensed phase to return to the gas phase
    Real k_r_{ 0 };
  };
}  // namespace micm
