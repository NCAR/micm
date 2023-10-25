// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#define TIMED_METHOD(assigned_increment, time_it, method, ...)                                 \
  {                                                                                            \
    if constexpr (time_it)                                                                     \
    {                                                                                          \
      auto start = std::chrono::high_resolution_clock::now();                                  \
      method(__VA_ARGS__);                                                                     \
      auto end = std::chrono::high_resolution_clock::now();                                    \
      assigned_increment += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start); \
    }                                                                                          \
    else                                                                                       \
    {                                                                                          \
      method(__VA_ARGS__);                                                                     \
    }                                                                                          \
  }

namespace micm
{
  // annonymous namespace to hide jacobian builder
  namespace {
    template<template<class> class SparseMatrixPolicy>
    SparseMatrixPolicy<double> build_jacobian(
      std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements,
      size_t number_of_grid_cells,
      size_t state_size
    )
    {
      auto builder = SparseMatrixPolicy<double>::create(state_size).number_of_blocks(number_of_grid_cells);
      for (auto& elem : nonzero_jacobian_elements)
        builder = builder.with_element(elem.first, elem.second);
      // Always include diagonal elements
      for (std::size_t i = 0; i < state_size; ++i)
        builder = builder.with_element(i, i);
      
      return SparseMatrixPolicy<double>(builder);
    }
  }
  //
  // RosenbrockSolverParameters
  //
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

  //
  // RosenbrockSolver
  //
  inline void SolverStats::Reset()
  {
    function_calls = 0;
    jacobian_updates = 0;
    number_of_steps = 0;
    accepted = 0;
    rejected = 0;
    decompositions = 0;
    solves = 0;
    singular = 0;
    total_forcing_time = std::chrono::nanoseconds::zero();
    total_jacobian_time = std::chrono::nanoseconds::zero();
    total_linear_factor_time = std::chrono::nanoseconds::zero();
    total_linear_solve_time = std::chrono::nanoseconds::zero();
  }

  inline std::string StateToString(const SolverState& state)
  {
    switch (state)
    {
      case SolverState::NotYetCalled: return "Not Yet Called";
      case SolverState::Running: return "Running";
      case SolverState::Converged: return "Converged";
      case SolverState::ConvergenceExceededMaxSteps: return "Convergence Exceeded Max Steps";
      case SolverState::StepSizeTooSmall: return "Step Size Too Small";
      case SolverState::RepeatedlySingularMatrix: return "Repeatedly Singular Matrix";
      case SolverState::NaNDetected: return "NaNDetected";
      default: return "Unknown";
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::RosenbrockSolver()
      : processes_(),
        parameters_(RosenbrockSolverParameters::three_stage_rosenbrock_parameters()),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters)
      : RosenbrockSolver(
            system,
            processes,
            parameters,
            [](const SparseMatrixPolicy<double>& matrix, double initial_value) -> LinearSolverPolicy {
              return LinearSolverPolicy{ matrix, initial_value };
            })
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters,
      const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver)
      : processes_(processes),
        parameters_(parameters),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
    std::map<std::string, std::size_t> variable_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;

    std::size_t index = 0;
    for (auto& name : system.UniqueNames())
      variable_map[name] = index++;

    // generate a state-vector reordering function to reduce fill-in in linear solver
    if (parameters_.reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_process_set = ProcessSet(processes, variable_map);
      auto unsorted_jac_elements = unsorted_process_set.NonZeroJacobianElements();
      MatrixPolicy<int> unsorted_jac_non_zeros(system.StateSize(), system.StateSize(), 0);
      for (auto& elem : unsorted_jac_elements)
        unsorted_jac_non_zeros[elem.first][elem.second] = 1;
      auto reorder_map = DiagonalMarkowitzReorder<MatrixPolicy>(unsorted_jac_non_zeros);
      state_reordering = [=](const std::vector<std::string>& variables, const std::size_t i)
      { return variables[reorder_map[i]]; };

      variable_map.clear();
      std::size_t index = 0;
      for (auto& name : system.UniqueNames(state_reordering))
        variable_map[name] = index++;
    }
    
    // setup the state_parameters
    std::vector<std::string> param_labels{};
    for (const auto& process : processes_)
      if (process.rate_constant_)
        for (auto& label : process.rate_constant_->CustomParameters())
          param_labels.push_back(label);

    process_set_ = ProcessSet(processes, variable_map);

    auto jacobian = build_jacobian<SparseMatrixPolicy>(
      process_set_.NonZeroJacobianElements(),
      parameters_.number_of_grid_cells_,
      system.StateSize()
    );

    std::vector<std::size_t> jacobian_diagonal_elements;
    for (std::size_t i = 0; i < jacobian[0].size(); ++i)
      jacobian_diagonal_elements.push_back(jacobian.VectorIndex(0, i, i));

    state_parameters_ = {
      .number_of_grid_cells_ = parameters_.number_of_grid_cells_,
      .number_of_rate_constants_ = processes_.size(),
      .variable_names_ = system.UniqueNames(state_reordering),
      .custom_rate_parameter_labels_ = param_labels,
      .jacobian_diagonal_elements_ = jacobian_diagonal_elements,
      .state_size_ = system.StateSize()
    };

    process_set_.SetJacobianFlatIds(jacobian);
    linear_solver_ = std::move(create_linear_solver(jacobian, 1.0e-30));
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy> RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::GetState() const
  {
    auto state = State<MatrixPolicy, SparseMatrixPolicy>{ state_parameters_ };
      
    state.jacobian_ = build_jacobian<SparseMatrixPolicy>(
      process_set_.NonZeroJacobianElements(),
      state_parameters_.number_of_grid_cells_,
      state_parameters_.state_size_
    );

    auto lu = linear_solver_.GetLUMatrices(state.jacobian_, 1.0e-30);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    state.lower_matrix_ = lower_matrix;
    state.upper_matrix_ = upper_matrix;

    return state;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  template<bool time_it>
  inline typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::SolverResult
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::Solve(
      double time_step,
      State<MatrixPolicy, SparseMatrixPolicy>& state) noexcept
  {
    typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::SolverResult result{};
    result.state_ = SolverState::Running;
    MatrixPolicy<double> Y(state.variables_);
    MatrixPolicy<double> Ynew(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> initial_forcing(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> forcing(Y.size(), Y[0].size(), 0.0);
    MatrixPolicy<double> temp(Y.size(), Y[0].size(), 0.0);
    std::vector<MatrixPolicy<double>> K{};
    const double h_max = parameters_.h_max_ == 0.0 ? time_step : std::min(time_step, parameters_.h_max_);
    const double h_start =
        parameters_.h_start_ == 0.0 ? std::max(parameters_.h_min_, delta_min_) : std::min(h_max, parameters_.h_start_);

    SolverStats stats;
    UpdateState(state);

    for (std::size_t i = 0; i < parameters_.stages_; ++i)
      K.push_back(MatrixPolicy<double>(Y.size(), Y[0].size(), 0.0));

    double present_time = 0.0;
    double H = std::min(std::max(std::abs(parameters_.h_min_), std::abs(h_start)), std::abs(h_max));

    if (std::abs(H) <= 10 * parameters_.round_off_)
      H = delta_min_;

    // TODO: the logic above this point should be moved to the constructor and should return an error
    //       if the parameters are invalid (e.g., h_min > h_max or h_start > h_max)

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_step + parameters_.round_off_) <= 0 && (result.state_ == SolverState::Running))
    {
      if (stats.number_of_steps > parameters_.max_number_of_steps_)
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

      // compute the forcing at the beginning of the current time
      TIMED_METHOD(stats.total_forcing_time, time_it, CalculateForcing, state.rate_constants_, Y, initial_forcing);
      stats.function_calls += 1;

      // compute the jacobian at the beginning of the current time
      TIMED_METHOD(stats.total_jacobian_time, time_it, CalculateJacobian, state.rate_constants_, Y, state.jacobian_);
      stats.jacobian_updates += 1;

      bool accepted = false;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        bool is_singular{ false };
        // Form and factor the rosenbrock ode jacobian
        TIMED_METHOD(stats.total_linear_factor_time, time_it, LinearFactor, H, parameters_.gamma_[0], is_singular, Y, stats, state);
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
              TIMED_METHOD(stats.total_forcing_time, time_it, CalculateForcing, state.rate_constants_, Ynew, forcing);
              stats.function_calls += 1;
            }
          }
          K[stage].AsVector().assign(forcing.AsVector().begin(), forcing.AsVector().end());
          for (uint64_t j = 0; j < stage; ++j)
          {
            auto HC = parameters_.c_[stage_combinations + j] / H;
            K[stage].ForEach([&](double& iKstage, double& iKj) { iKstage += HC * iKj; }, K[j]);
          }
          temp.AsVector().assign(K[stage].AsVector().begin(), K[stage].AsVector().end());
          TIMED_METHOD(stats.total_linear_solve_time, time_it, linear_solver_.template Solve<MatrixPolicy>, temp, K[stage], state.lower_matrix_, state.upper_matrix_);
          stats.solves += 1;
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

        stats.number_of_steps += 1;

        // Check the error magnitude and adjust step size
        if (std::isnan(error))
        {
          Y.AsVector().assign(Ynew.AsVector().begin(), Ynew.AsVector().end());
          result.state_ = SolverState::NaNDetected;
          break;
        }
        else if ((error < 1) || (H < parameters_.h_min_))
        {
          stats.accepted += 1;
          present_time = present_time + H;
          Y.AsVector().assign(Ynew.AsVector().begin(), Ynew.AsVector().end());
          Hnew = std::max(parameters_.h_min_, std::min(Hnew, h_max));
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
          if (stats.accepted >= 1)
          {
            stats.rejected += 1;
          }
        }
      }
    }

    if (result.state_ == SolverState::Running)
    {
      result.state_ = SolverState::Converged;
    }

    result.final_time_ = present_time;
    result.stats_ = stats;
    result.result_ = Y;
    return result;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing)
  {
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
    process_set_.AddForcingTerms<MatrixPolicy>(rate_constants, number_densities, forcing);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const
    requires(!VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    for (auto& elem : jacobian.AsVector())
      elem = -elem;
    for (std::size_t i_block = 0; i_block < jacobian.size(); ++i_block)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_block * jacobian.FlatBlockSize());
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::AlphaMinusJacobian(
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
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::CalculateJacobian(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      SparseMatrixPolicy<double>& jacobian)
  {
    std::fill(jacobian.AsVector().begin(), jacobian.AsVector().end(), 0.0);
    process_set_.AddJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(rate_constants, number_densities, jacobian);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::UpdateState(State<MatrixPolicy, SparseMatrixPolicy>& state)
  {
    Process::UpdateState(processes_, state);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::LinearFactor(
      double& H,
      const double gamma,
      bool& singular,
      const MatrixPolicy<double>& number_densities,
      SolverStats& stats, 
      State<MatrixPolicy, SparseMatrixPolicy>& state)
  {
    auto jacobian = state.jacobian_;
    uint64_t n_consecutive = 0;
    singular = false;
    while (true)
    {
      double alpha = 1 / (H * gamma);
      AlphaMinusJacobian(jacobian, alpha);
      if (parameters_.check_singularity_)
      {
        linear_solver_.Factor(jacobian, state.lower_matrix_, state.upper_matrix_, singular);
      } else {
        singular = false;
        linear_solver_.Factor(jacobian, state.lower_matrix_, state.upper_matrix_);
      }
      singular = false;  // TODO This should be evaluated in some way
      stats.decompositions += 1;
      if (!singular)
        break;
      stats.singular += 1;
      if (++n_consecutive > 5)
        break;
      H /= 2;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
  inline double RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>::NormalizedError(
      const MatrixPolicy<double>& Y,
      const MatrixPolicy<double>& Ynew,
      const MatrixPolicy<double>& errors)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    auto _y = Y.AsVector();
    auto _ynew = Ynew.AsVector();
    auto _errors = errors.AsVector();
    size_t N = Y.AsVector().size();

    double error = 0;

    for (size_t i = 0; i < N; ++i)
    {
      double ymax = std::max(std::abs(_y[i]), std::abs(_ynew[i]));
      double scale = parameters_.absolute_tolerance_ + parameters_.relative_tolerance_ * ymax;
      error += std::pow(_errors[i] / scale, 2);
    }

    double error_min_ = 1.0e-10;
    return std::max(std::sqrt(error / N), error_min_);
  }

}  // namespace micm
