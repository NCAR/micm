// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <vector>

namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy>
  class AbstractBackwardEuler;

  /// @brief BackwardEuler solver parameters
  struct BackwardEulerSolverParameters
  {
    template<class RatesPolicy, class LinearSolverPolicy>
    using SolverType = AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy>;

    double small_{ 1.0e-40 };
    double h_start_{ 0.0 };
    size_t max_number_of_steps_{ 11 };
    // The time step reductions are used to determine the time step after a failed solve
    // This default set of time step reductions is used by CAM-Chem
    std::array<double, 5> time_step_reductions_{ 0.5, 0.5, 0.5, 0.5, 0.1 };
  };

}  // namespace micm