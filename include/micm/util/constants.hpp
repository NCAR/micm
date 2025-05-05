// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace micm
{
  namespace constants
  {
    static constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;                      // J K^{-1}
    static constexpr double AVOGADRO_CONSTANT = 6.02214076e23;                      // # mol^{-1}
    static constexpr double GAS_CONSTANT = BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT;  // J K^{-1} mol^{-1}
  }  // namespace constants
}  // namespace micm
