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

    void print() const
    {
      std::cout << "stages_: " << stages_ << std::endl;
      std::cout << "upper_limit_tolerance_: " << upper_limit_tolerance_ << std::endl;
      std::cout << "max_number_of_steps_: " << max_number_of_steps_ << std::endl;
      std::cout << "round_off_: " << round_off_ << std::endl;
      std::cout << "factor_min_: " << factor_min_ << std::endl;
      std::cout << "factor_max_: " << factor_max_ << std::endl;
      std::cout << "rejection_factor_decrease_: " << rejection_factor_decrease_ << std::endl;
      std::cout << "safety_factor_: " << safety_factor_ << std::endl;
      std::cout << "h_min_: " << h_min_ << std::endl;
      std::cout << "h_max_: " << h_max_ << std::endl;
      std::cout << "h_start_: " << h_start_ << std::endl;
      std::cout << "new_function_evaluation_: ";
      for (bool val : new_function_evaluation_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "estimator_of_local_order_: " << estimator_of_local_order_ << std::endl;
      std::cout << "a_: ";
      for (double val : a_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "c_: ";
      for (double val : c_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "m_: ";
      for (double val : m_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "e_: ";
      for (double val : e_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "alpha_: ";
      for (double val : alpha_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "gamma_: ";
      for (double val : gamma_)
        std::cout << val << " ";
      std::cout << std::endl;
      std::cout << "absolute_tolerance_: " << absolute_tolerance_ << std::endl;
      std::cout << "relative_tolerance_: " << relative_tolerance_ << std::endl;
      std::cout << "number_of_grid_cells_: " << number_of_grid_cells_ << std::endl;
    }

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

  RosenbrockSolverParameters RosenbrockSolverParameters::two_stage_rosenbrock_parameters(
      size_t number_of_grid_cells,
      bool reorder_state)
  {
    // an L-stable method, 2 stages, order 2

    RosenbrockSolverParameters parameters;
    double g = 1.0 + 1.0 / std::sqrt(2.0);

    parameters.stages_ = 2;

    parameters.a_.fill(0);
    parameters.a_[0] = 1.0 / g;

    parameters.c_.fill(0);
    parameters.c_[0] = -2.0 / g;

    // Both stages require a new function evaluation
    parameters.new_function_evaluation_.fill(true);

    parameters.m_.fill(0);
    parameters.m_[0] = (3.0) / (2.0 * g);
    parameters.m_[1] = (1.0) / (2.0 * g);

    parameters.e_.fill(0);
    parameters.e_[0] = 1.0 / (2.0 * g);
    parameters.e_[1] = 1.0 / (2.0 * g);

    parameters.estimator_of_local_order_ = 2.0;

    parameters.alpha_.fill(0);
    parameters.alpha_[0] = 0.0;
    parameters.alpha_[1] = 1.0;

    parameters.gamma_.fill(0);
    parameters.gamma_[0] = g;
    parameters.gamma_[1] = -g;

    parameters.number_of_grid_cells_ = number_of_grid_cells;
    parameters.reorder_state_ = reorder_state;

    return parameters;
  }

  RosenbrockSolverParameters RosenbrockSolverParameters::three_stage_rosenbrock_parameters(
      size_t number_of_grid_cells,
      bool reorder_state)
  {
    // an L-stable method, 3 stages, order 3, 2 function evaluations
    //
    // original formaulation for three stages:
    // Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997.
    // Benchmarking stiff ode solvers for atmospheric chemistry problems II: Rosenbrock solvers.
    // Atmospheric Environment 31, 3459–3472. https://doi.org/10.1016/S1352-2310(97)83212-8
    RosenbrockSolverParameters parameters;

    parameters.stages_ = 3;

    parameters.a_.fill(0);
    parameters.a_[0] = 1;
    parameters.a_[1] = 1;
    parameters.a_[2] = 0;

    parameters.c_.fill(0);
    parameters.c_[0] = -0.10156171083877702091975600115545e+01;
    parameters.c_[1] = 0.40759956452537699824805835358067e+01;
    parameters.c_[2] = 0.92076794298330791242156818474003e+01;

    parameters.new_function_evaluation_.fill(false);
    parameters.new_function_evaluation_[0] = true;
    parameters.new_function_evaluation_[1] = true;
    parameters.new_function_evaluation_[2] = false;
    parameters.m_.fill(0);
    parameters.m_[0] = 0.1e+01;
    parameters.m_[1] = 0.61697947043828245592553615689730e+01;
    parameters.m_[2] = -0.42772256543218573326238373806514;

    parameters.e_.fill(0);
    parameters.e_[0] = 0.5;
    parameters.e_[1] = -0.29079558716805469821718236208017e+01;
    parameters.e_[2] = 0.22354069897811569627360909276199;

    parameters.estimator_of_local_order_ = 3;

    parameters.alpha_.fill(0);
    parameters.alpha_[0] = 0;
    parameters.alpha_[1] = 0.43586652150845899941601945119356;
    parameters.alpha_[2] = 0.43586652150845899941601945119356;

    parameters.gamma_.fill(0);
    parameters.gamma_[0] = 0.43586652150845899941601945119356;
    parameters.gamma_[1] = 0.24291996454816804366592249683314;
    parameters.gamma_[2] = 0.21851380027664058511513169485832e+01;

    parameters.number_of_grid_cells_ = number_of_grid_cells;
    parameters.reorder_state_ = reorder_state;

    return parameters;
  }

  RosenbrockSolverParameters RosenbrockSolverParameters::four_stage_rosenbrock_parameters(
      size_t number_of_grid_cells,
      bool reorder_state)
  {
    // L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
    // L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
    //
    //  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
    //  EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
    //  SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
    //  SPRINGER-VERLAG (1990)
    RosenbrockSolverParameters parameters;

    parameters.stages_ = 4;

    parameters.a_.fill(0);
    parameters.a_[0] = 0.2000000000000000E+01;
    parameters.a_[1] = 0.1867943637803922E+01;
    parameters.a_[2] = 0.2344449711399156;
    parameters.a_[3] = parameters.a_[1];
    parameters.a_[4] = parameters.a_[2];
    parameters.a_[5] = 0.0;

    parameters.c_.fill(0);
    parameters.c_[0] = -0.7137615036412310E+01;
    parameters.c_[1] = 0.2580708087951457E+01;
    parameters.c_[2] = 0.6515950076447975;
    parameters.c_[3] = -0.2137148994382534E+01;
    parameters.c_[4] = -0.3214669691237626;
    parameters.c_[5] = -0.6949742501781779;

    parameters.new_function_evaluation_.fill(false);
    parameters.new_function_evaluation_[0] = true;
    parameters.new_function_evaluation_[1] = true;
    parameters.new_function_evaluation_[2] = true;
    parameters.new_function_evaluation_[3] = false;

    parameters.m_.fill(0);
    parameters.m_[0] = 0.2255570073418735E+01;
    parameters.m_[1] = 0.2870493262186792;
    parameters.m_[2] = 0.4353179431840180;
    parameters.m_[3] = 0.1093502252409163E+01;

    parameters.e_.fill(0);
    parameters.e_[0] = -0.2815431932141155;
    parameters.e_[1] = -0.7276199124938920E-01;
    parameters.e_[2] = -0.1082196201495311;
    parameters.e_[3] = -0.1093502252409163E+01;

    parameters.estimator_of_local_order_ = 4.0;

    parameters.alpha_.fill(0);
    parameters.alpha_[0] = 0.0;
    parameters.alpha_[1] = 0.1145640000000000E+01;
    parameters.alpha_[2] = 0.6552168638155900;
    parameters.alpha_[3] = parameters.alpha_[2];

    parameters.gamma_.fill(0);
    parameters.gamma_[0] = 0.5728200000000000;
    parameters.gamma_[1] = -0.1769193891319233E+01;
    parameters.gamma_[2] = 0.7592633437920482;
    parameters.gamma_[3] = -0.1049021087100450;

    parameters.number_of_grid_cells_ = number_of_grid_cells;
    parameters.reorder_state_ = reorder_state;

    return parameters;
  }

  RosenbrockSolverParameters RosenbrockSolverParameters::four_stage_differential_algebraic_rosenbrock_parameters(
      size_t number_of_grid_cells,
      bool reorder_state)
  {
    // A STIFFLY-STABLE METHOD, 4 stages, order 3
    RosenbrockSolverParameters parameters;

    // Set the number of stages
    parameters.stages_ = 4;

    parameters.a_.fill(0.0);
    parameters.a_[0] = 0.0;
    parameters.a_[1] = 2.0;
    parameters.a_[2] = 0.0;
    parameters.a_[3] = 2.0;
    parameters.a_[4] = 0.0;
    parameters.a_[5] = 1.0;

    parameters.c_.fill(0.0);
    parameters.c_[0] = 4.0;
    parameters.c_[1] = 1.0;
    parameters.c_[2] = -1.0;
    parameters.c_[3] = 1.0;
    parameters.c_[4] = -1.0;
    parameters.c_[5] = -(8.0 / 3.0);

    parameters.new_function_evaluation_.fill(false);
    parameters.new_function_evaluation_[0] = true;
    parameters.new_function_evaluation_[2] = true;
    parameters.new_function_evaluation_[3] = true;

    parameters.m_.fill(0.0);
    parameters.m_[0] = 2.0;
    parameters.m_[2] = 1.0;
    parameters.m_[3] = 1.0;

    parameters.e_.fill(0.0);
    parameters.e_[3] = 1.0;

    parameters.estimator_of_local_order_ = 3.0;

    parameters.alpha_.fill(0.0);
    parameters.alpha_[2] = 1.0;
    parameters.alpha_[3] = 1.0;

    parameters.gamma_.fill(0.0);
    parameters.gamma_[0] = 0.5;
    parameters.gamma_[1] = 1.5;

    parameters.number_of_grid_cells_ = number_of_grid_cells;
    parameters.reorder_state_ = reorder_state;

    return parameters;
  }

  RosenbrockSolverParameters RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters(
      size_t number_of_grid_cells,
      bool reorder_state)
  {
    // STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
    //
    // E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
    // EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
    // SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
    // SPRINGER-VERLAG (1996)

    RosenbrockSolverParameters parameters;

    parameters.stages_ = 6;

    parameters.alpha_.fill(0.0);
    parameters.alpha_[0] = 0.000;
    parameters.alpha_[1] = 0.386;
    parameters.alpha_[2] = 0.210;
    parameters.alpha_[3] = 0.630;
    parameters.alpha_[4] = 1.000;
    parameters.alpha_[5] = 1.000;

    parameters.gamma_.fill(0.0);
    parameters.gamma_[0] = 0.2500000000000000;
    parameters.gamma_[1] = -0.1043000000000000;
    parameters.gamma_[2] = 0.1035000000000000;
    parameters.gamma_[3] = -0.3620000000000023E-01;
    parameters.gamma_[4] = 0.0;
    parameters.gamma_[5] = 0.0;

    parameters.a_.fill(0.0);
    parameters.a_[0] = 0.1544000000000000E+01;
    parameters.a_[1] = 0.9466785280815826;
    parameters.a_[2] = 0.2557011698983284;
    parameters.a_[3] = 0.3314825187068521E+01;
    parameters.a_[4] = 0.2896124015972201E+01;
    parameters.a_[5] = 0.9986419139977817;
    parameters.a_[6] = 0.1221224509226641E+01;
    parameters.a_[7] = 0.6019134481288629E+01;
    parameters.a_[8] = 0.1253708332932087E+02;
    parameters.a_[9] = -0.6878860361058950;
    parameters.a_[10] = parameters.a_[6];
    parameters.a_[11] = parameters.a_[7];
    parameters.a_[12] = parameters.a_[8];
    parameters.a_[13] = parameters.a_[9];
    parameters.a_[14] = 1.0;

    parameters.c_.fill(0.0);
    parameters.c_[0] = -0.5668800000000000E+01;
    parameters.c_[1] = -0.2430093356833875E+01;
    parameters.c_[2] = -0.2063599157091915;
    parameters.c_[3] = -0.1073529058151375;
    parameters.c_[4] = -0.9594562251023355E+01;
    parameters.c_[5] = -0.2047028614809616E+02;
    parameters.c_[6] = 0.7496443313967647E+01;
    parameters.c_[7] = -0.1024680431464352E+02;
    parameters.c_[8] = -0.3399990352819905E+02;
    parameters.c_[9] = 0.1170890893206160E+02;
    parameters.c_[10] = 0.8083246795921522E+01;
    parameters.c_[11] = -0.7981132988064893E+01;
    parameters.c_[12] = -0.3152159432874371E+02;
    parameters.c_[13] = 0.1631930543123136E+02;
    parameters.c_[14] = -0.6058818238834054E+01;

    parameters.m_.fill(0.0);
    parameters.m_[0] = parameters.a_[6];
    parameters.m_[1] = parameters.a_[7];
    parameters.m_[2] = parameters.a_[8];
    parameters.m_[3] = parameters.a_[9];
    parameters.m_[4] = 1.0;
    parameters.m_[5] = 1.0;

    parameters.e_.fill(0.0);
    parameters.e_[5] = 1.0;

    parameters.new_function_evaluation_.fill(true);

    parameters.estimator_of_local_order_ = 4.0;

    parameters.number_of_grid_cells_ = number_of_grid_cells;
    parameters.reorder_state_ = reorder_state;

    return parameters;
  }

  /// @brief An implementation of the Chapman mechnanism solver
  ///
  /// The template parameter is the type of matrix to use
  template<template<class> class MatrixPolicy = Matrix, template<class> class SparseMatrixPolicy = SparseMatrix>
  class RosenbrockSolver
  {
   public:
    enum class SolverState
    {
      NotYetCalled,
      Converged,
      ConvergenceExceededMaxSteps,
      StepSizeTooSmall,
      RepeatedlySingularMatrix,
      NaNDetected
    };

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
      std::string State(const RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState& state) const;
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

    const System system_;
    const std::vector<Process> processes_;
    RosenbrockSolverParameters parameters_;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering_;
    ProcessSet process_set_;
    SolverStats stats_;
    SparseMatrixPolicy<double> jacobian_;
    LinearSolver<double, SparseMatrixPolicy> linear_solver_;
    std::vector<std::size_t> jacobian_diagonal_elements_;
    size_t N_{};

    static constexpr double delta_min_ = 1.0e-5;

    /// @brief Default constructor
    RosenbrockSolver();

    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    RosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters);

    virtual ~RosenbrockSolver();

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
    void CalculateForcing(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& number_densities,
        MatrixPolicy<double>& forcing);

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    inline void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
      requires(!VectorizableSparse<SparseMatrixPolicy<double>>);
    inline void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
      requires(VectorizableSparse<SparseMatrixPolicy<double>>);

    /// @brief Update the rate constants for the environment state
    /// @param state The current state of the chemical system
    void UpdateState(State<MatrixPolicy>& state);

    /// @brief Compute the derivative of the forcing w.r.t. each chemical, the jacobian
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param jacobian The matrix of partial derivatives
    /// @return The jacobian
    void CalculateJacobian(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& number_densities,
        SparseMatrixPolicy<double>& jacobian);

    /// @brief Prepare the linear solver
    /// @param H time step (seconds)
    /// @param gamma time step factor for specific rosenbrock method
    /// @param singular indicates if the matrix is singular
    /// @param number_densities constituent concentration (molec/cm^3)
    /// @param rate_constants Rate constants for each process (molecule/cm3)^(n-1) s-1
    void LinearFactor(
        double& H,
        const double gamma,
        bool& singular,
        const MatrixPolicy<double>& number_densities,
        const MatrixPolicy<double>& rate_constants);

   protected:
    /// @brief Computes the scaled norm of the vector errors
    /// @param original_number_densities the original number densities
    /// @param new_number_densities the new number densities
    /// @param errors The computed errors
    /// @return
    double NormalizedError(
        MatrixPolicy<double> original_number_densities,
        MatrixPolicy<double> new_number_densities,
        MatrixPolicy<double> errors);
  };

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverStats::Reset()
  {
    function_calls = 0;
    jacobian_updates = 0;
    number_of_steps = 0;
    accepted = 0;
    rejected = 0;
    decompositions = 0;
    solves = 0;
    singular = 0;
    total_steps = 0;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  std::string RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverStats::State(
      const RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState& state) const
  {
    switch (state)
    {
      case RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState::NotYetCalled: return "Not Yet Called";
      case RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState::Converged: return "Converged";
      case RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState::ConvergenceExceededMaxSteps:
        return "Convergence Exceeded Max Steps";
      case RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState::StepSizeTooSmall: return "Step Size Too Small";
      case RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverState::RepeatedlySingularMatrix:
        return "Repeatedly Singular Matrix";
      default: return "Unknown";
    }
    return "";
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::RosenbrockSolver()
      : system_(),
        processes_(),
        parameters_(RosenbrockSolverParameters::three_stage_rosenbrock_parameters()),
        process_set_(),
        stats_(),
        jacobian_(),
        linear_solver_(),
        jacobian_diagonal_elements_(),
        N_(system_.StateSize() * parameters_.number_of_grid_cells_)
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters)
      : system_(system),
        processes_(processes),
        parameters_(parameters),
        state_reordering_(),
        process_set_(),
        stats_(),
        jacobian_(),
        linear_solver_(),
        jacobian_diagonal_elements_(),
        N_(system_.StateSize() * parameters_.number_of_grid_cells_)
  {
    // generate a state-vector reordering function to reduce fill-in in linear solver
    if (parameters_.reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_process_set = ProcessSet(processes, GetState());
      auto unsorted_jac_elements = unsorted_process_set.NonZeroJacobianElements();
      MatrixPolicy<int> unsorted_jac_non_zeros(system_.StateSize(), system_.StateSize(), 0);
      for (auto& elem : unsorted_jac_elements)
        unsorted_jac_non_zeros[elem.first][elem.second] = 1;
      auto reorder_map = DiagonalMarkowitzReorder<MatrixPolicy>(unsorted_jac_non_zeros);
      state_reordering_ = [=](const std::vector<std::string>& variables, const std::size_t i)
      { return variables[reorder_map[i]]; };
    }
    process_set_ = ProcessSet(processes, GetState());
    auto builder =
        SparseMatrixPolicy<double>::create(system_.StateSize()).number_of_blocks(parameters_.number_of_grid_cells_);
    auto jac_elements = process_set_.NonZeroJacobianElements();
    for (auto& elem : jac_elements)
      builder = builder.with_element(elem.first, elem.second);
    // Always include diagonal elements
    for (std::size_t i = 0; i < system_.StateSize(); ++i)
      builder = builder.with_element(i, i);

    jacobian_ = builder;
    linear_solver_ = LinearSolver<double, SparseMatrixPolicy>(jacobian_, 1.0e-30);
    process_set_.SetJacobianFlatIds(jacobian_);
    for (std::size_t i = 0; i < jacobian_[0].size(); ++i)
      jacobian_diagonal_elements_.push_back(jacobian_.VectorIndex(0, i, i));
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::~RosenbrockSolver()
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline State<MatrixPolicy> RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::GetState() const
  {
    std::vector<std::string> param_labels{};
    for (const auto& process : processes_)
      if (process.rate_constant_)
        for (auto& label : process.rate_constant_->CustomParameters())
          param_labels.push_back(label);
    return State<MatrixPolicy>{ micm::StateParameters{ .state_variable_names_ = system_.UniqueNames(state_reordering_),
                                                       .custom_rate_parameter_labels_ = param_labels,
                                                       .number_of_grid_cells_ = parameters_.number_of_grid_cells_,
                                                       .number_of_rate_constants_ = processes_.size() } };
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverResult
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::Solve(double time_step, State<MatrixPolicy>& state) noexcept
  {
    typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::SolverResult result{};
    result.state_ = SolverState::Converged;
    MatrixPolicy<double> Y(state.variables_);
    MatrixPolicy<double> Ynew(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> initial_forcing(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> forcing(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> temp(Y.size(), Y[0].size(), 0.0);
    std::vector<MatrixPolicy<double>> K{};

    parameters_.h_max_ = time_step;
    parameters_.h_start_ = std::max(parameters_.h_min_, delta_min_);

    stats_.Reset();
    UpdateState(state);

    for (std::size_t i = 0; i < parameters_.stages_; ++i)
      K.push_back(MatrixPolicy<double>(Y.size(), Y[0].size(), 0.0));

    double present_time = 0.0;
    double H =
        std::min(std::max(std::abs(parameters_.h_min_), std::abs(parameters_.h_start_)), std::abs(parameters_.h_max_));

    if (std::abs(H) <= 10 * parameters_.round_off_)
      H = delta_min_;

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_step + parameters_.round_off_) <= 0 && (result.state_ == SolverState::Converged))
    {
      if (stats_.number_of_steps > parameters_.max_number_of_steps_)
      {
        result.state_ = SolverState::ConvergenceExceededMaxSteps;
        break;
      }

      if (((present_time + 0.1 * H) == present_time) || (H <= parameters_.round_off_))
      {
        result.state_ = SolverState::StepSizeTooSmall;
        break;
      }

      //  Limit H if necessary to avoid going beyond the specified chemistry time step
      H = std::min(H, std::abs(time_step - present_time));

      // compute the concentrations at the current time
      CalculateForcing(state.rate_constants_, Y, initial_forcing);

      // compute the jacobian at the current time
      CalculateJacobian(state.rate_constants_, Y, jacobian_);

      bool accepted = false;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        bool is_singular{ false };
        // Form and factor the rosenbrock ode jacobian
        LinearFactor(H, parameters_.gamma_[0], is_singular, Y, state.rate_constants_);
        if (is_singular)
        {
          result.state_ = SolverState::RepeatedlySingularMatrix;
          break;
        }

        // Compute the stages
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
        {
          double stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
          if (stage == 0)
          {
            forcing = initial_forcing;
          }
          else
          {
            if (parameters_.new_function_evaluation_[stage])
            {
              Ynew.AsVector().assign(Y.AsVector().begin(), Y.AsVector().end());
              for (uint64_t j = 0; j < stage; ++j)
              {
                auto a = parameters_.a_[stage_combinations + j];
                Ynew.ForEach([&](double& iYnew, double& iKj) { iYnew += a * iKj; }, K[j]);
              }
              CalculateForcing(state.rate_constants_, Ynew, forcing);
            }
          }
          K[stage].AsVector().assign(forcing.AsVector().begin(), forcing.AsVector().end());
          for (uint64_t j = 0; j < stage; ++j)
          {
            auto HC = parameters_.c_[stage_combinations + j] / H;
            K[stage].ForEach([&](double& iKstage, double& iKj) { iKstage += HC * iKj; }, K[j]);
          }
          temp.AsVector().assign(K[stage].AsVector().begin(), K[stage].AsVector().end());
          linear_solver_.template Solve<MatrixPolicy>(temp, K[stage]);
          stats_.solves += 1;
        }

        // Compute the new solution
        Ynew.AsVector().assign(Y.AsVector().begin(), Y.AsVector().end());
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Ynew.ForEach([&](double& iYnew, double& iKstage) { iYnew += parameters_.m_[stage] * iKstage; }, K[stage]);

        // Compute the error estimation
        MatrixPolicy<double> Yerror(Y.size(), Y[0].size(), 0);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Yerror.ForEach([&](double& iYerror, double& iKstage) { iYerror += parameters_.e_[stage] * iKstage; }, K[stage]);

        auto error = NormalizedError(Y, Ynew, Yerror);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
        double fac = std::min(
                              parameters_.factor_max_,
                              std::max(
                                  parameters_.factor_min_,
                                  parameters_.safety_factor_ / std::pow(error, 1 / parameters_.estimator_of_local_order_)));
        double Hnew = H * fac;

        // Check the error magnitude and adjust step size
        stats_.number_of_steps += 1;
        stats_.total_steps += 1;

        if (std::isnan(error)) {
          result.state_ = SolverState::NaNDetected;
          break;
        }
        else if ((error < 1) || (H < parameters_.h_min_))
        {
          stats_.accepted += 1;
          present_time = present_time + H;
          Y.AsVector().assign(Ynew.AsVector().begin(), Ynew.AsVector().end());
          Hnew = std::max(parameters_.h_min_, std::min(Hnew, parameters_.h_max_));
          if (reject_last_h)
          {
            // No step size increase after a rejected step
            Hnew = std::min(Hnew, H);
          }
          reject_last_h = false;
          reject_more_h = false;
          H = Hnew;
          accepted = true;
        }
        else
        {
          // Reject step
          if (reject_more_h)
          {
            Hnew = H * parameters_.rejection_factor_decrease_;
          }
          reject_more_h = reject_last_h;
          reject_last_h = true;
          H = Hnew;
          if (stats_.accepted >= 1)
          {
            stats_.rejected += 1;
          }
        }
      }
    }

    result.final_time_ = present_time;
    result.stats_ = stats_;
    result.result_ = std::move(Y);
    return result;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing)
  {
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
    process_set_.AddForcingTerms<MatrixPolicy>(rate_constants, number_densities, forcing);
    stats_.function_calls += 1;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const
    requires(!VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    for (auto& elem : jacobian.AsVector())
      elem = -elem;
    for (std::size_t i_block = 0; i_block < jacobian.size(); ++i_block)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_block * jacobian.FlatBlockSize());
      for (const auto& i_elem : jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const
    requires(VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    const std::size_t n_cells = jacobian.GroupVectorSize();
    for (auto& elem : jacobian.AsVector())
      elem = -elem;
    for (std::size_t i_group = 0; i_group < jacobian.NumberOfGroups(jacobian.size()); ++i_group)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_group * jacobian.GroupSize(jacobian.FlatBlockSize()));
      for (const auto& i_elem : jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::CalculateJacobian(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      SparseMatrixPolicy<double>& jacobian)
  {
    std::fill(jacobian.AsVector().begin(), jacobian.AsVector().end(), 0.0);
    process_set_.AddJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(rate_constants, number_densities, jacobian);
    stats_.jacobian_updates += 1;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::UpdateState(State<MatrixPolicy>& state)
  {
    Process::UpdateState(processes_, state);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::LinearFactor(
      double& H,
      const double gamma,
      bool& singular,
      const MatrixPolicy<double>& number_densities,
      const MatrixPolicy<double>& rate_constants)
  {
    // TODO: invesitage this function. The fortran equivalent appears to have a bug.
    // From my understanding the fortran do loop would only ever do one iteration and is equivalent to what's below
    SparseMatrixPolicy<double> jacobian = jacobian_;
    uint64_t n_consecutive = 0;
    singular = true;
    while (true)
    {
      double alpha = 1 / (H * gamma);
      AlphaMinusJacobian(jacobian, alpha);
      linear_solver_.Factor(jacobian);
      singular = false;  // TODO This should be evaluated in some way
      stats_.decompositions += 1;
      if (!singular)
        break;
      stats_.singular += 1;
      if (++n_consecutive > 5)
        break;
      H /= 2;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline double RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>::NormalizedError(
      MatrixPolicy<double> Y,
      MatrixPolicy<double> Ynew,
      MatrixPolicy<double> errors)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7
    MatrixPolicy<double> maxs(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> scale(Y.size(), Y[0].size(), 0.0);

    maxs.ForEach([&](double& imax, double& iY, double& iYnew) { imax = std::max(std::abs(iY), std::abs(iYnew)); }, Y, Ynew);

    scale.ForEach(
        [&](double& iscale, double& imax)
        { iscale = parameters_.absolute_tolerance_ + parameters_.relative_tolerance_ * imax; },
        maxs);

    double sum = 0;
    errors.ForEach([&](double& ierror, double& iscale) { sum += std::pow(ierror / iscale, 2); }, scale);

    double error_min_ = 1.0e-10;
    return std::max(std::sqrt(sum / N_), error_min_);
  }
}  // namespace micm