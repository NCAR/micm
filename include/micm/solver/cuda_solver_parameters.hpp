/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include "rosenbrock_solver_parameters.hpp"

#include <micm/solver/cuda_rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

namespace micm
{
  /// @brief Parameters for the CUDA Rosenbrock solver
  struct CudaRosenbrockSolverParameters : public RosenbrockSolverParameters
  {
    template<class RatesPolicy, class LinearSolverPolicy>
    using SolverType = CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>;

    /// @brief Constructor from base class
    /// @param base
    CudaRosenbrockSolverParameters(const RosenbrockSolverParameters& base)
        : RosenbrockSolverParameters(base)
    {
    }
  };
}  // namespace micm