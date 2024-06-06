// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "rosenbrock_solver_parameters.hpp"
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/jit_rosenbrock.hpp>

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
}