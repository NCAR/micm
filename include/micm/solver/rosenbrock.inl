// Copyright (C) 2023-2024 National Center for Atmospheric Research
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
    total_update_state_time = std::chrono::nanoseconds::zero();
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

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::RosenbrockSolver()
      : processes_(),
        parameters_(RosenbrockSolverParameters::three_stage_rosenbrock_parameters()),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters)
      : RosenbrockSolver(
            system,
            processes,
            parameters,
            [](const SparseMatrixPolicy<double>& matrix, double initial_value) -> LinearSolverPolicy {
              return LinearSolverPolicy{ matrix, initial_value };
            },
            [](const std::vector<Process>& processes, const std::map<std::string, std::size_t>& variable_map) -> ProcessSetPolicy {
              return ProcessSetPolicy{ processes, variable_map };
            })
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters,
      const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
      const std::function<ProcessSetPolicy(const std::vector<Process>&, const std::map<std::string, std::size_t>&)> create_process_set)
      : processes_(processes),
        parameters_(parameters),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
    std::map<std::string, std::size_t> variable_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;

    std::size_t index = 0;
    auto used_species = ProcessSetPolicy::SpeciesUsed(processes);
    auto available_species = system.UniqueNames();
    std::sort(available_species.begin(), available_species.end());
    std::set<std::string> unused_species;
    std::set_difference(available_species.begin(), available_species.end(), used_species.begin(), used_species.end(), std::inserter(unused_species, unused_species.begin()));
    if (unused_species.size() > 0 && !parameters_.ignore_unused_species_)
    {
      std::string err_msg = "Unused species in chemical system:";
      for (auto& species: unused_species)
        err_msg += " '" + species + "'";
      err_msg += ". Set solver parameter ignore_unused_species_ to allow unused species in solve.";
      throw std::runtime_error(err_msg);
    }
    for (auto& name : system.UniqueNames())
      variable_map[name] = index++;

    // generate a state-vector reordering function to reduce fill-in in linear solver
    if (parameters_.reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_process_set = std::move(create_process_set(processes, variable_map));
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

    process_set_ = std::move(create_process_set(processes, variable_map));

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
      .nonzero_jacobian_elements_ = process_set_.NonZeroJacobianElements()
    };

    process_set_.SetJacobianFlatIds(jacobian);
    linear_solver_ = std::move(create_linear_solver(jacobian, 1.0e-30));
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy> RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::GetState() const
  {
    return State<MatrixPolicy, SparseMatrixPolicy>{ state_parameters_ };   
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  template<bool time_it>
  inline typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::SolverResult
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::Solve(
      double time_step,
      State<MatrixPolicy, SparseMatrixPolicy>& state) noexcept
  {
    typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::SolverResult result{};
    result.state_ = SolverState::Running;
    // reset the upper, lower matrix. Repeated calls without zeroing these matrices can lead to lack of convergence
    auto& lower = state.lower_matrix_.AsVector();
    auto& upper = state.upper_matrix_.AsVector();
    for(auto& l : lower) l = 0;
    for(auto& u : upper) u = 0;
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
    TIMED_METHOD(stats.total_update_state_time, time_it, UpdateState, state);

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

      // compute the negative jacobian at the beginning of the current time
      TIMED_METHOD(stats.total_jacobian_time, time_it, CalculateNegativeJacobian, state.rate_constants_, Y, state.jacobian_);
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
                Ynew.ForEach([&](double& iYnew, const double& iKj) { iYnew += a * iKj; }, K[j]);
              }
              TIMED_METHOD(stats.total_forcing_time, time_it, CalculateForcing, state.rate_constants_, Ynew, forcing);
              stats.function_calls += 1;
            }
          }
          K[stage].AsVector().assign(forcing.AsVector().begin(), forcing.AsVector().end());
          for (uint64_t j = 0; j < stage; ++j)
          {
            auto HC = parameters_.c_[stage_combinations + j] / H;
            K[stage].ForEach([&](double& iKstage, const double& iKj) { iKstage += HC * iKj; }, K[j]);
          }
          temp.AsVector().assign(K[stage].AsVector().begin(), K[stage].AsVector().end());
          TIMED_METHOD(stats.total_linear_solve_time, time_it, linear_solver_.template Solve<MatrixPolicy>, temp, K[stage], state.lower_matrix_, state.upper_matrix_);
          stats.solves += 1;
        }

        // Compute the new solution
        Ynew.AsVector().assign(Y.AsVector().begin(), Y.AsVector().end());
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Ynew.ForEach([&](double& iYnew, const double& iKstage) { iYnew += parameters_.m_[stage] * iKstage; }, K[stage]);

        // Compute the error estimation
        MatrixPolicy<double> Yerror(Y.size(), Y[0].size(), 0);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Yerror.ForEach([&](double& iYerror, const double& iKstage) { iYerror += parameters_.e_[stage] * iKstage; }, K[stage]);

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

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing)
  {
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
    process_set_.template AddForcingTerms<MatrixPolicy>(rate_constants, number_densities, forcing);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const
    requires(!VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    for (std::size_t i_block = 0; i_block < jacobian.size(); ++i_block)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_block * jacobian.FlatBlockSize());
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const
    requires(VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    const std::size_t n_cells = jacobian.GroupVectorSize();
     
    for (std::size_t i_group = 0; i_group < jacobian.NumberOfGroups(jacobian.size()); ++i_group)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_group * jacobian.GroupSize(jacobian.FlatBlockSize()));
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::CalculateNegativeJacobian(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      SparseMatrixPolicy<double>& jacobian)
  {
    std::fill(jacobian.AsVector().begin(), jacobian.AsVector().end(), 0.0);
    process_set_.template SubtractJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(rate_constants, number_densities, jacobian);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::UpdateState(State<MatrixPolicy, SparseMatrixPolicy>& state)
  {
    Process::UpdateState(processes_, state);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::LinearFactor(
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

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
  inline double RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::NormalizedError(
      const MatrixPolicy<double>& Y,
      const MatrixPolicy<double>& Ynew,
      const MatrixPolicy<double>& errors) const
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
