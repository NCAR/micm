// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/backward_euler_dae_solver_parameters.hpp>
#include <micm/solver/backward_euler_dae_temporary_variables.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
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

  /// @brief An implementation of the fully implicit backward euler method with algebraic constraints
  /// 
  /// This solver extends the standard backward Euler method to handle Differential Algebraic Equations (DAEs)
  /// which consist of both differential equations and algebraic constraints.
  /// 
  /// The DAE system has the form:
  /// - Differential equations: dy/dt = f(t, y)
  /// - Algebraic constraints: 0 = g(t, y)
  /// 
  /// Usage example:
  /// @code
  /// auto params = micm::BackwardEulerDAESolverParameters();
  /// 
  /// // Add a mass conservation constraint: A + B + C = constant
  /// params.algebraic_constraints_.push_back([initial_mass](double time, 
  ///     const std::vector<double>& variables, const std::vector<double>& rate_constants) {
  ///   if (variables.size() >= 3) {
  ///     double current_mass = variables[0] + variables[1] + variables[2]; // A + B + C
  ///     return current_mass - initial_mass; // Should be zero for conservation
  ///   }
  ///   return 0.0;
  /// });
  /// 
  /// // Add a time-dependent constraint: species A = 2*t
  /// params.algebraic_constraints_.push_back([](double time, 
  ///     const std::vector<double>& variables, const std::vector<double>& rate_constants) {
  ///   if (variables.size() >= 1) {
  ///     return variables[0] - 2.0 * time; // A should equal 2*time
  ///   }
  ///   return 0.0;
  /// });
  /// @endcode
  /// 
  template<class RatesPolicy, class LinearSolverPolicy>
  class AbstractBackwardEulerDAE
  {
   public:
    LinearSolverPolicy linear_solver_;
    RatesPolicy rates_;

    /// @brief Solver parameters typename
    using ParametersType = BackwardEulerDAESolverParameters;

    /// @brief Default constructor
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    AbstractBackwardEulerDAE(
        LinearSolverPolicy&& linear_solver,
        RatesPolicy&& rates,
        auto& jacobian,
        const size_t number_of_species)
        : linear_solver_(std::move(linear_solver)),
          rates_(std::move(rates))
    {
    }

    AbstractBackwardEulerDAE(const AbstractBackwardEulerDAE&) = delete;
    AbstractBackwardEulerDAE& operator=(const AbstractBackwardEulerDAE&) = delete;
    AbstractBackwardEulerDAE(AbstractBackwardEulerDAE&&) = default;
    AbstractBackwardEulerDAE& operator=(AbstractBackwardEulerDAE&&) = default;

    virtual ~AbstractBackwardEulerDAE() = default;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @param state The state to advance
    /// @return result of the solver (success or failure, and statistics)
    SolverResult Solve(double time_step, auto& state, const BackwardEulerDAESolverParameters& parameters) const;

    /// @brief Determines whether the residual is small enough to stop the
    ///        internal solver iteration
    /// @param residual The residual to check
    /// @param state The current state being solved for
    /// @return true if the residual is small enough to stop the iteration
    template<class DenseMatrixPolicy>
    static bool IsConverged(
        const BackwardEulerDAESolverParameters& parameters,
        const DenseMatrixPolicy& residual,
        const DenseMatrixPolicy& Yn1,
        const std::vector<double>& absolute_tolerance,
        double relative_tolerance)
      requires(!VectorizableDense<DenseMatrixPolicy>);
    template<class DenseMatrixPolicy>
    static bool IsConverged(
        const BackwardEulerDAESolverParameters& parameters,
        const DenseMatrixPolicy& residual,
        const DenseMatrixPolicy& Yn1,
        const std::vector<double>& absolute_tolerance,
        double relative_tolerance)
      requires(VectorizableDense<DenseMatrixPolicy>);

   private:
    /// @brief Evaluate algebraic constraints
    /// @param time Current time
    /// @param state Current state
    /// @param parameters Solver parameters containing constraint functions
    /// @param constraint_values Output vector for constraint values
    template<class DenseMatrixPolicy>
    void EvaluateConstraints(
        double time,
        const auto& state,
        const BackwardEulerDAESolverParameters& parameters,
        std::vector<double>& constraint_values) const;

   public:
    /// @brief Check if algebraic constraints are satisfied
    /// @param constraint_values Vector of constraint values
    /// @param tolerance Tolerance for constraint satisfaction
    /// @return true if all constraints are satisfied within tolerance
    bool AreConstraintsSatisfied(
        const std::vector<double>& constraint_values,
        double tolerance = 1e-12) const;
  };

}  // namespace micm

#include "backward_euler_dae.inl"