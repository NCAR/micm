/* Copyright (C) 2023 National Center for Atmospheric Research,
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
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <string>
#include <vector>

namespace micm
{

  /// @brief Rosenbrock solver parameters
  struct RosenbrockSolverParameters
  {
    size_t stages_{};
    size_t upper_limit_tolerance_{};
    size_t max_number_of_steps_{ 1000 };

    double round_off_{ std::numeric_limits<double>::epsilon() };  // Unit roundoff (1+round_off)>1
    double factor_min_{ 0.2 };                                    // solver step size minimum boundary
    double factor_max_{ 6 };                                      // solver step size maximum boundary
    double rejection_factor_decrease_{ 0.1 };                     // used to decrease the step after 2 successive rejections
    double safety_factor_{ 0.9 };                                 // safety factor in new step size computation

    double h_min_{ 0 };        // step size min
    double h_max_{ 0.5 };      // step size max
    double h_start_{ 0.005 };  // step size start

    // Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
    // or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    std::array<bool, 6>
        new_function_evaluation_{};  // which steps reuse the previous iterations evaluation or do a new evaluation

    double estimator_of_local_order_{};  // the minumum between the main and the embedded scheme orders plus one

    //  The coefficient matrices A and C are strictly lower triangular.
    //  The lower triangular (subdiagonal) elements are stored in row-wise order:
    //  A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
    //  The general mapping formula is:
    //      A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
    //      C(i,j) = ros_C( (i-1)*(i-2)/2 + j )
    std::array<double, 15> a_{};  // coefficient matrix a
    std::array<double, 15> c_{};  // coefficient matrix c
    std::array<double, 6> m_{};   // coefficients for new step evaluation
    std::array<double, 6> e_{};   // error estimation coefficients

    // Y_stage_i ~ Y( T + H*Alpha_i )
    std::array<double, 6> alpha_{};
    // Gamma_i = \sum_j  gamma_{i,j}
    std::array<double, 6> gamma_{};

    double absolute_tolerance_{ 1e-3 };
    double relative_tolerance_{ 1e-4 };

    size_t number_of_grid_cells_{ 1 };  // Number of grid cells to solve simultaneously
    bool reorder_state_{ true };        // Reorder state during solver construction to minimize LU fill-in

    // Print RosenbrockSolverParameters to console
    void print() const;

    static RosenbrockSolverParameters two_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    static RosenbrockSolverParameters three_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    static RosenbrockSolverParameters four_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);

    static RosenbrockSolverParameters four_stage_differential_algebraic_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    static RosenbrockSolverParameters six_stage_differential_algebraic_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);

   private:
    RosenbrockSolverParameters() = default;
  };

  /// @brief The final state the solver was in after the Solve function finishes
  enum class SolverState
  {
    NotYetCalled,
    Running,
    Converged,
    ConvergenceExceededMaxSteps,
    StepSizeTooSmall,
    RepeatedlySingularMatrix,
    NaNDetected
  };

  std::string StateToString(const SolverState& state);

  /// @brief An implementation of the Rosenbrock ODE solver
  ///
  /// The template parameter is the type of matrix to use
  template<template<class> class MatrixPolicy = Matrix, template<class> class SparseMatrixPolicy = SparseMatrix>
  class RosenbrockSolver
  {
   public:
    struct SolverStats
    {
      uint64_t function_calls{};    // Nfun
      uint64_t jacobian_updates{};  // Njac
      uint64_t number_of_steps{};   // Nstp
      uint64_t accepted{};          // Nacc
      uint64_t rejected{};          // Nrej
      uint64_t decompositions{};    // Ndec
      uint64_t solves{};            // Nsol
      uint64_t singular{};          // Nsng
      uint64_t total_steps{};       // Ntotstp

      void Reset();
    };

    struct [[nodiscard]] SolverResult
    {
      /// @brief The new state computed by the solver
      MatrixPolicy<double> result_{};
      /// @brief The finals state the solver was in
      SolverState state_ = SolverState::NotYetCalled;
      /// @brief A collection of runtime state for this call of the solver
      SolverStats stats_{};
      /// @brief The final time the solver iterated to
      double final_time_{};
    };

    System system_;
    std::vector<Process> processes_;
    RosenbrockSolverParameters parameters_;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering_;
    ProcessSet process_set_;
    SolverStats stats_;
    SparseMatrixPolicy<double> jacobian_;
    LinearSolver<double, SparseMatrixPolicy> linear_solver_;
    std::vector<std::size_t> jacobian_diagonal_elements_;
    size_t N_{};

    static constexpr double delta_min_ = 1.0e-6;

    /// @brief Default constructor
    RosenbrockSolver();

    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    RosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters);

    virtual ~RosenbrockSolver() = default;

    /// @brief Returns a state object for use with the solver
    /// @return A object that can hold the full state of the chemical system
    State<MatrixPolicy> GetState() const;

    /// @brief Advances the given step over the specified time step
    /// @param time_step Time [s] to advance the state by
    /// @return A struct containing results and a status code
    SolverResult Solve(double time_step, State<MatrixPolicy>& state) noexcept;

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
    void UpdateState(State<MatrixPolicy>& state);

    /// @brief Compute the derivative of the forcing w.r.t. each chemical, the jacobian
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param jacobian The matrix of partial derivatives
    virtual void CalculateJacobian(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& number_densities,
        SparseMatrixPolicy<double>& jacobian);

    /// @brief Prepare the linear solver
    /// @param H time step (seconds)
    /// @param gamma time step factor for specific rosenbrock method
    /// @param singular indicates if the matrix is singular
    /// @param number_densities constituent concentration (molec/cm^3)
    /// @param rate_constants Rate constants for each process (molecule/cm3)^(n-1) s-1
    void LinearFactor(double& H, const double gamma, bool& singular, const MatrixPolicy<double>& number_densities);

   protected:
    /// @brief Computes the scaled norm of the vector errors
    /// @param Y the original vector
    /// @param Ynew the new vector
    /// @param errors The computed errors
    /// @return
    double
    NormalizedError(const MatrixPolicy<double>& y, const MatrixPolicy<double>& y_new, const MatrixPolicy<double>& errors);
  };

}  // namespace micm

#include "rosenbrock.inl"
