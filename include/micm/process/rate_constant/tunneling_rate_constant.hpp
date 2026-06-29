// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace micm
{
  struct TunnelingRateConstantParameters
  {
    /// @brief Pre-exponential factor [(mol m−3)^(−(𝑛−1)) s−1]
    double A_ = 1.0;
    /// @brief Linear temperature-dependent parameter [K]
    double B_ = 0.0;
    /// @brief Cubed temperature-dependent parameter [K^3]
    double C_ = 0.0;
  };
}  // namespace micm
