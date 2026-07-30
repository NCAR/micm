// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <array>
#include <cstddef>
#include <iostream>
#include <vector>

namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
  class AbstractBackwardEuler;

  /// @brief Backward Euler solver parameters
  struct BackwardEulerSolverParameters
  {
    template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
    using SolverType = AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>;

    Real small_{ 1.0e-40 };
    Real h_start_{ 0.0 };
    Index max_number_of_steps_{ 11 };
    // The time step reductions are used to determine the time step after a failed solve
    // This default set of time step reductions is used by CAM-Chem
    std::array<Real, 5> time_step_reductions_{ 0.5, 0.5, 0.5, 0.5, 0.1 };
  };

}  // namespace micm