// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <vector>

namespace micm
{

  /// @brief BackwardEuler solver parameters
  struct BackwardEulerSolverParameters
  {
    std::vector<double> absolute_tolerance_;
    double small{ 1.0e-40 };
    size_t max_number_of_steps_{ 11 };
    std::array<double, 5> time_step_reductions{ 0.5, 0.5, 0.5, 0.5, 0.1 };
  };

}  // namespace micm