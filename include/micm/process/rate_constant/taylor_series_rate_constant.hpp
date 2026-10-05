// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <cstddef>

namespace micm
{
  struct TaylorSeriesRateConstantParameters
  {
    /// @brief Maximum number of Taylor series coefficients supported
    static constexpr Index MAX_COEFFICIENTS = 16;
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
    /// @brief Taylor coefficients for the series expansion.
    ///        Only the first n_coefficients_ entries are used.
    Real coefficients_[MAX_COEFFICIENTS] = { 1.0 };
    /// @brief Number of active Taylor coefficients [1, MAX_COEFFICIENTS]
    Index n_coefficients_{ 1 };
  };
}  // namespace micm
