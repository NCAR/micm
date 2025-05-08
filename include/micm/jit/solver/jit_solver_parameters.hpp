/* Copyright (C) 2023-2025 University Corporation for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/jit/solver/jit_rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

namespace micm
{
  /// @brief Parameters for the JIT Rosenbrock solver
  struct JitRosenbrockSolverParameters : public RosenbrockSolverParameters
  {
    template<class RatesPolicy, class LinearSolverPolicy>
    using SolverType = JitRosenbrockSolver<RatesPolicy, LinearSolverPolicy>;

    /// @brief Constructor from base class
    /// @param base
    JitRosenbrockSolverParameters(const RosenbrockSolverParameters& base)
        : RosenbrockSolverParameters(base)
    {
    }
  };
}  // namespace micm