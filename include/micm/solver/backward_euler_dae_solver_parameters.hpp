// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy>
  class AbstractBackwardEulerDAE;

  /// @brief BackwardEuler DAE solver parameters with algebraic constraints
  struct BackwardEulerDAESolverParameters
  {
    template<class RatesPolicy, class LinearSolverPolicy>
    using SolverType = AbstractBackwardEulerDAE<RatesPolicy, LinearSolverPolicy>;

    double small_{ 1.0e-40 };
    double h_start_{ 0.0 };
    size_t max_number_of_steps_{ 11 };
    // The time step reductions are used to determine the time step after a failed solve
    // This default set of time step reductions is used by CAM-Chem
    std::array<double, 5> time_step_reductions_{ 0.5, 0.5, 0.5, 0.5, 0.1 };
    
    /// @brief Algebraic constraint functions
    /// Each function takes (time, state_variables, rate_constants) and returns the constraint value
    /// The constraint should equal zero when satisfied
    /// Using a type-erased approach for constraints
    std::vector<std::function<double(double, const std::vector<double>&, const std::vector<double>&)>> algebraic_constraints_;
  };

}  // namespace micm