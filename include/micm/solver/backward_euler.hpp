/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/state.hpp>
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
  class BackwardEuler
  {
   public:
    /// @brief Default constructor
    BackwardEuler();

    /// @brief Builds a backward eulder solver for the given system and processes
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    BackwardEuler(const System& system, const std::vector<Process>& processes);

    virtual ~BackwardEuler() = default;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @param state The state to advance
    /// @param linear_solver The linear solver to use
    /// @param process_set The process set to use
    /// @param processes The collection of chemical processes that will be applied during solving
    /// @param jacobian_diagonal_elements The diagonal elements of the jacobian matrix
    /// @return A struct containing results and a status code
    void Solve(
        double time_step,
        auto& state,
        auto linear_solver,
        auto process_set,
        const std::vector<micm::Process>& processes,
        auto jacobian_diagonal_elements);
  };

}  // namespace micm

#include "backward_euler.inl"
