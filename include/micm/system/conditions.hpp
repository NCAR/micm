// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/constants.hpp>
#include <micm/util/types.hpp>

namespace micm
{
  /// @brief Environemental conditions
  struct Conditions
  {
    Real temperature_{ 0.0 };  // K
    Real pressure_{ 0.0 };     // Pa
    Real air_density_{ 0.0 };  // mol m-3

    void CalculateIdealAirDensity();
  };

  inline void Conditions::CalculateIdealAirDensity()
  {
    air_density_ = pressure_ / (constants::GAS_CONSTANT * temperature_);
  }
}  // namespace micm
