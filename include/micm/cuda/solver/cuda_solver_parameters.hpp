/* Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/cuda/solver/cuda_rosenbrock.hpp>
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