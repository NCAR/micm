/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

namespace micm
{
  /// @brief Environemental conditions
  struct Conditions
  {
    double temperature_{ 0.0 };  // K
    double pressure_{ 0.0 };     // Pa;
    double air_density_{ 1.0 };  // mol m-3
  };
}  // namespace micm
