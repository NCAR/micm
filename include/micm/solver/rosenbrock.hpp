// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Much of this solver was formulated and implemented from this book:
// Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
// edition. ed. Springer, Berlin ; New York. The source code for many (all?) of the solvers in that book can be found here:
// http://www.unige.ch/~hairer/software.html
//
// Some extensions to the rosenbrock solver formulated there were formulated in this paper
// Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. Benchmarking stiff ode solvers for
// atmospheric chemistry problems II: Rosenbrock solvers. Atmospheric Environment 31, 3459–3472.
// https://doi.org/10.1016/S1352-2310(97)83212-8
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/rosenbrock_temporary_variables.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <algorithm>
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

  /// @brief An implementation of the Rosenbrock ODE solver
  /// @tparam RatesPolicy Calculator of forcing and Jacobian terms
  /// @tparam LinearSolverPolicy Linear solver
  /// @tparam Derived Implementation of the Rosenbock solver
  ///
  /// This implements the Curiously Recurring Template Pattern to allow
  /// the AlphaMinusJacobian and NormalizedError functions to be implemented
  /// in extending classes and called from the base class Solve() function.
  /// https://en.cppreference.com/w/cpp/language/crtp
  ///
  template<class RatesPolicy, class LinearSolverPolicy, class Dervied>
  class AbstractRosenbrockSolver
  {
   public:
    LinearSolverPolicy linear_solver_;
    RatesPolicy rates_;

    static constexpr double DELTA_MIN = 1.0e-6;

    /// @brief Solver parameters typename
    using ParametersType = RosenbrockSolverParameters;

    /// @brief Default constructor
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    /// Note: This constructor is not intended to be used directly. Instead, use the SolverBuilder to create a solver
    AbstractRosenbrockSolver(
        LinearSolverPolicy&& linear_solver,
        RatesPolicy&& rates,
        auto& jacobian,
        const size_t number_of_species)
        : linear_solver_(std::move(linear_solver)),
          rates_(std::move(rates))
    {
    }

    AbstractRosenbrockSolver(const AbstractRosenbrockSolver&) = delete;
    AbstractRosenbrockSolver& operator=(const AbstractRosenbrockSolver&) = delete;
    AbstractRosenbrockSolver(AbstractRosenbrockSolver&&) = default;
    AbstractRosenbrockSolver& operator=(AbstractRosenbrockSolver&&) = default;

    virtual ~AbstractRosenbrockSolver() = default;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @return A struct containing results and a status code
    SolverResult Solve(double time_step, auto& state, const RosenbrockSolverParameters& parameters) const noexcept;

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    template<class SparseMatrixPolicy>
    void AlphaMinusJacobian(
        SparseMatrixPolicy& jacobian,
        std::vector<std::size_t>& jacobian_diagonal_elements,
        const double& alpha) const
      requires(!VectorizableSparse<SparseMatrixPolicy>);
    template<class SparseMatrixPolicy>
    void AlphaMinusJacobian(
        SparseMatrixPolicy& jacobian,
        std::vector<std::size_t>& jacobian_diagonal_elements,
        const double& alpha) const
      requires(VectorizableSparse<SparseMatrixPolicy>);

    /// @brief Perform the LU decomposition of the matrix
    /// @param alpha The alpha value
    /// @param number_densities The number densities
    /// @param stats The solver stats
    /// @param state The state
    void LinearFactor(const double alpha, SolverStats& stats, auto& state) const;

    /// @brief Computes the scaled norm of the vector errors
    /// @param y the original vector
    /// @param y_new the new vector
    /// @param errors The computed errors
    /// @return
    template<class DenseMatrixPolicy>
    double NormalizedError(
        const DenseMatrixPolicy& y,
        const DenseMatrixPolicy& y_new,
        const DenseMatrixPolicy& errors,
        auto& state) const
      requires(!VectorizableDense<DenseMatrixPolicy>);
    template<class DenseMatrixPolicy>
    double NormalizedError(
        const DenseMatrixPolicy& y,
        const DenseMatrixPolicy& y_new,
        const DenseMatrixPolicy& errors,
        auto& state) const
      requires(VectorizableDense<DenseMatrixPolicy>);
  };  // end of Abstract Rosenbrock Solver

  template<class RatesPolicy, class LinearSolverPolicy>
  class RosenbrockSolver
      : public AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, RosenbrockSolver<RatesPolicy, LinearSolverPolicy>>
  {
   public:
    /// @brief Default constructor
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    /// @param jacobian Jacobian matrix
    ///
    /// Note: This constructor is not intended to be used directly. Instead, use the SolverBuilder to create a solver
    RosenbrockSolver(LinearSolverPolicy&& linear_solver, RatesPolicy&& rates, auto& jacobian, const size_t number_of_species)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, RosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(
              std::move(linear_solver),
              std::move(rates),
              jacobian,
              number_of_species)
    {
    }
    RosenbrockSolver(const RosenbrockSolver&) = delete;
    RosenbrockSolver& operator=(const RosenbrockSolver&) = delete;
    RosenbrockSolver(RosenbrockSolver&&) = default;
    RosenbrockSolver& operator=(RosenbrockSolver&&) = default;
  };

}  // namespace micm

#include "rosenbrock.inl"
