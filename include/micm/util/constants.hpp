// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

namespace micm::constants
{
  static constexpr Real BOLTZMANN_CONSTANT = 1.380649e-23;                      // J K^{-1}
  static constexpr Real AVOGADRO_CONSTANT = 6.02214076e23;                      // # mol^{-1}
  static constexpr Real GAS_CONSTANT = BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT;  // J K^{-1} mol^{-1}
}  // namespace micm::constants
