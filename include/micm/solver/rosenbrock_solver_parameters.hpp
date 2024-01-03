// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <iostream>

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

    double h_min_{ 0.0 };  // step size min [s]
    double h_max_{
      0.0
    };  // step size max [s] (if zero or greater than the solver time-step, the time-step passed to the solver will be used)
    double h_start_{ 0.0 };  // step size start [s] (if zero, 1.0e-6 will be used, if greater than h_max, h_max will be used)

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

    size_t number_of_grid_cells_{ 1 };     // Number of grid cells to solve simultaneously
    bool reorder_state_{ true };           // Reorder state during solver construction to minimize LU fill-in
    bool check_singularity_{ false };      // Check for singular A matrix in linear solve of A x = b
    bool ignore_unused_species_{ false };  // Allow unused species to be included in state and solve

    // Print RosenbrockSolverParameters to console
    void print() const;

    /// @brief an L-stable method, 2 stages, order 2
    /// @param number_of_grid_cells
    /// @param reorder_state
    /// @return
    static RosenbrockSolverParameters two_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    /// @brief an L-stable method, 3 stages, order 3, 2 function evaluations
    /// @param number_of_grid_cells
    /// @param reorder_state
    /// @return
    static RosenbrockSolverParameters three_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    /// @brief L-stable rosenbrock method of order 4, with 4 stages
    /// @param number_of_grid_cells
    /// @param reorder_state
    /// @return
    static RosenbrockSolverParameters four_stage_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    /// @brief A stiffly-stable method, 4 stages, order 3
    /// @param number_of_grid_cells
    /// @param reorder_state
    /// @return
    static RosenbrockSolverParameters four_stage_differential_algebraic_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);
    /// @brief stiffly-stable rosenbrock method of order 4, with 6 stages
    /// @param number_of_grid_cells
    /// @param reorder_state
    /// @return
    static RosenbrockSolverParameters six_stage_differential_algebraic_rosenbrock_parameters(
        size_t number_of_grid_cells = 1,
        bool reorder_state = true);

   private:
    RosenbrockSolverParameters() = default;
  };

  inline void RosenbrockSolverParameters::print() const
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

  inline RosenbrockSolverParameters RosenbrockSolverParameters::two_stage_rosenbrock_parameters(
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

  inline RosenbrockSolverParameters RosenbrockSolverParameters::three_stage_rosenbrock_parameters(
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

  inline RosenbrockSolverParameters RosenbrockSolverParameters::four_stage_rosenbrock_parameters(
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

  inline RosenbrockSolverParameters RosenbrockSolverParameters::four_stage_differential_algebraic_rosenbrock_parameters(
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

  inline RosenbrockSolverParameters RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters(
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

}  // namespace micm