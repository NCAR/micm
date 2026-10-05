// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

namespace micm
{
  struct ArrheniusRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    Real A_{ 1 };
    /// @brief Unitless exponential factor
    Real B_{ 0 };
    /// @brief Activation threshold, expected to be the negative activation energy divided by the boltzman constant
    ///        [-E_a / k_b), K]
    Real C_{ 0 };
    /// @brief A factor that determines temperature dependence [K]
    Real D_{ 300 };
    /// @brief A factor that determines pressure dependence [Pa-1]
    Real E_{ 0 };
  };
}  // namespace micm
