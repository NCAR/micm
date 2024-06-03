/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/state.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace micm
{

  /// @brief An implementation of the fully implicit backward euler method
  template<class LinearSolverPolicy, class ProcessSetPolicy>
  class BackwardEuler
  {
    BackwardEulerSolverParameters parameters_;
    LinearSolverPolicy linear_solver_;
    ProcessSetPolicy process_set_;
    std::vector<std::size_t> jacobian_diagonal_elements_;
    std::vector<micm::Process> processes_;

   public:
    /// @brief Solver parameters typename
    using ParametersType = BackwardEulerSolverParameters;

    /// @brief Default constructor
    BackwardEuler(
        BackwardEulerSolverParameters parameters,
        LinearSolverPolicy linear_solver,
        ProcessSetPolicy process_set,
        std::vector<std::size_t> jacobian_diagonal_elements,
        std::vector<micm::Process>& processes)
        : parameters_(parameters),
          linear_solver_(linear_solver),
          process_set_(process_set),
          jacobian_diagonal_elements_(jacobian_diagonal_elements),
          processes_(processes)
    {
    }

    virtual ~BackwardEuler() = default;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @param state The state to advance
    /// @return result of the solver (success or failure, and statistics)
    SolverResult Solve(double time_step, auto& state);
  };

}  // namespace micm

#include "backward_euler.inl"
