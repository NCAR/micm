/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Much of this solver was formulated and implemented from this book:
 * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
 * edition. ed. Springer, Berlin ; New York. The source code for many (all?) of the solvers in that book can be found here:
 * http://www.unige.ch/~hairer/software.html
 *
 * Some extensions to the rosenbrock solver formulated there were formulated in this paper
 * Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. Benchmarking stiff ode solvers for
 * atmospheric chemistry problems II: Rosenbrock solvers. Atmospheric Environment 31, 3459–3472.
 * https://doi.org/10.1016/S1352-2310(97)83212-8
 *
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <string>
#include <vector>

namespace micm
{

  /// @brief The final state the solver was in after the Solve function finishes
  enum class SolverState
  {
    /// @brief This is the initial value at the start of the Solve function
    NotYetCalled,
    /// @brief This is only used for control flow in the Solve function
    Running,
    /// @brief A successful integration will have this value
    Converged,
    /// @brief If the number of steps exceeds the maximum value on the solver parameter, this value will be returned
    ConvergenceExceededMaxSteps,
    /// @brief Very stiff systems will likely result in a step size that is not useable for the solver
    StepSizeTooSmall,
    /// @brief Matrices that are singular more than once will set this value. At present, this should never be returned
    RepeatedlySingularMatrix,
    /// @brief Mostly this value is returned by systems that tend toward chemical explosions
    NaNDetected
  };

  std::string StateToString(const SolverState& state);

  struct SolverStats
  {
    /// @brief The number of forcing function calls
    uint64_t function_calls{};
    /// @brief The number of jacobian function calls
    uint64_t jacobian_updates{};
    /// @brief The total number of internal time steps taken
    uint64_t number_of_steps{};
    /// @brief The number of accepted integrations
    uint64_t accepted{};
    /// @brief The number of rejected integrations
    uint64_t rejected{};
    /// @brief The number of LU decompositions
    uint64_t decompositions{};
    /// @brief The number of linear solves
    uint64_t solves{};
    /// @brief The number of times a singular matrix is detected. For now, this will always be zero as we assume the matrix
    /// is never singular
    uint64_t singular{};
    /// @brief The cumulative amount of time spent updating the state (including rate constant calculations)
    std::chrono::duration<double, std::nano> total_update_state_time{};
    /// @brief The cumulative amount of time spent calculating the forcing function
    std::chrono::duration<double, std::nano> total_forcing_time{};
    /// @brief The cumulative amount of time spent calculating the jacobian
    std::chrono::duration<double, std::nano> total_jacobian_time{};
    /// @brief The cumulative amount of time spent calculating the linear factorization
    std::chrono::duration<double, std::nano> total_linear_factor_time{};
    /// @brief The cumulative amount of time spent calculating the linear solve
    std::chrono::duration<double, std::nano> total_linear_solve_time{};

    /// @brief Set all member variables to zero
    void Reset();
  };

  /// @brief An implementation of the Rosenbrock ODE solver
  ///
  /// The template parameter is the type of matrix to use
  template<
      template<class> class MatrixPolicy = Matrix,
      template<class> class SparseMatrixPolicy = StandardSparseMatrix,
      class LinearSolverPolicy = LinearSolver<double, SparseMatrixPolicy>,
      class ProcessSetPolicy = ProcessSet>
  class RosenbrockSolver
  {
   public:
    struct [[nodiscard]] SolverResult
    {
      /// @brief The new state computed by the solver
      MatrixPolicy<double> result_{};
      /// @brief The final state the solver was in
      SolverState state_ = SolverState::NotYetCalled;
      /// @brief A collection of runtime state for this call of the solver
      SolverStats stats_{};
      /// @brief The final time the solver iterated to
      double final_time_{};
    };

    std::vector<Process> processes_;
    RosenbrockSolverParameters parameters_;
    StateParameters state_parameters_;
    ProcessSetPolicy process_set_;
    LinearSolverPolicy linear_solver_;

    static constexpr double delta_min_ = 1.0e-6;

    /// @brief Default constructor
    RosenbrockSolver();

    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    /// @param parameters Rosenbrock algorithm parameters
    RosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters);

    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters,
    ///        with a specific function provided to create the linear solver
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    /// @param parameters Rosenbrock algorithm parameters
    /// @param create_linear_solver Function that will be used to create a linear solver instance
    /// @param create_process_set Function that will be used to create a process set instance
    RosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters,
        const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
        const std::function<ProcessSetPolicy(const std::vector<Process>&, const std::map<std::string, std::size_t>&)>
            create_process_set);

    virtual ~RosenbrockSolver() = default;

    /// @brief Returns a state object for use with the solver
    /// @return A object that can hold the full state of the chemical system
    virtual State<MatrixPolicy, SparseMatrixPolicy> GetState() const;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @return A struct containing results and a status code
    template<bool time_it = false>
    SolverResult Solve(double time_step, State<MatrixPolicy, SparseMatrixPolicy>& state) noexcept;

    /// @brief Calculate a chemical forcing
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param forcing Vector of forcings for the current conditions
    virtual void CalculateForcing(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& number_densities,
        MatrixPolicy<double>& forcing);

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
        requires(!VectorizableSparse<SparseMatrixPolicy<double>>);
    void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
        requires(VectorizableSparse<SparseMatrixPolicy<double>>);

    /// @brief Update the rate constants for the environment state
    /// @param state The current state of the chemical system
    void UpdateState(State<MatrixPolicy, SparseMatrixPolicy>& state);

    /// @brief Compute the derivative of the forcing w.r.t. each chemical, and return the negative jacobian
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param jacobian The matrix of negative partial derivatives
    virtual void CalculateNegativeJacobian(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& number_densities,
        SparseMatrixPolicy<double>& jacobian);

    /// @brief Prepare the linear solver
    /// @param H time step (seconds)
    /// @param gamma time step factor for specific rosenbrock method
    /// @param singular indicates if the matrix is singular
    /// @param number_densities constituent concentration (molec/cm^3)
    /// @param state The current State
    void LinearFactor(
        double& H,
        const double gamma,
        bool& singular,
        const MatrixPolicy<double>& number_densities,
        SolverStats& stats,
        State<MatrixPolicy, SparseMatrixPolicy>& state);

   protected:
    /// @brief Computes the scaled norm of the vector errors
    /// @param y the original vector
    /// @param y_new the new vector
    /// @param errors The computed errors
    /// @return
    double NormalizedError(
        const MatrixPolicy<double>& y,
        const MatrixPolicy<double>& y_new,
        const MatrixPolicy<double>& errors) const;
  };

}  // namespace micm

#include "rosenbrock.inl"
