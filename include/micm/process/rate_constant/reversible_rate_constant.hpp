// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace micm
{
  struct ReversibleRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    double A_{ 1 };
    /// @brief Activation threshold [K]
    double C_{ 0 };
    /// @brief Reverse rate constant [s−1], indicating how fast the species
    ///        leaves the condensed phase to return to the gas phase
    double k_r_{ 0 };
  };
}  // namespace micm
