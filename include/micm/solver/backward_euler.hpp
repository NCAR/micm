// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/backward_euler_temporary_variables.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>

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
  template<class RatesPolicy, class LinearSolverPolicy>
  class BackwardEuler
  {
   public:
    BackwardEulerSolverParameters parameters_;
    LinearSolverPolicy linear_solver_;
    RatesPolicy rates_;
    std::vector<std::size_t> jacobian_diagonal_elements_;

    /// @brief Solver parameters typename
    using ParametersType = BackwardEulerSolverParameters;

    /// @brief Default constructor
    /// @param parameters Solver parameters
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    /// @param jacobian Jacobian matrix
    BackwardEuler(
        const BackwardEulerSolverParameters& parameters,
        LinearSolverPolicy&& linear_solver,
        RatesPolicy&& rates,
        auto& jacobian,
        const size_t number_of_species)
        : parameters_(parameters),
          linear_solver_(std::move(linear_solver)),
          rates_(std::move(rates)),
          jacobian_diagonal_elements_(jacobian.DiagonalIndices(0))
    {
    }

    BackwardEuler(const BackwardEuler&) = delete;
    BackwardEuler& operator=(const BackwardEuler&) = delete;
    BackwardEuler(BackwardEuler&&) = default;
    BackwardEuler& operator=(BackwardEuler&&) = default;

    virtual ~BackwardEuler() = default;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @param state The state to advance
    /// @return result of the solver (success or failure, and statistics)
    SolverResult Solve(double time_step, auto& state) const;

    /// @brief Determines whether the residual is small enough to stop the
    ///        internal solver iteration
    /// @param parameters Solver parameters
    /// @param residual The residual to check
    /// @param state The current state being solved for
    /// @return true if the residual is small enough to stop the iteration
    template<class DenseMatrixPolicy>
    static bool IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1, const std::vector<double>& absolute_tolerance, double relative_tolerance) requires(!VectorizableDense<DenseMatrixPolicy>);
    template<class DenseMatrixPolicy>
    static bool IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1, const std::vector<double>& absolute_tolerance, double relative_tolerance) requires(VectorizableDense<DenseMatrixPolicy>);
  };

}  // namespace micm

#include "backward_euler.inl"
