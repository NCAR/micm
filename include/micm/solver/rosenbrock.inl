// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

enum class MicmRosenbrockErrc
{
  UnusedSpecies = 1,  // Unused species present in the chemical system
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmRosenbrockErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class RosenbrockErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return "MICM Rosenbrock";
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmRosenbrockErrc>(ev))
      {
        case MicmRosenbrockErrc::UnusedSpecies:
          return "Unused species present in the chemical system. Use the ignore_unused_species_ parameter to allow unused "
                 "species in the solve.";
        default: return "Unknown error";
      }
    }
  };

  const RosenbrockErrorCategory rosenbrockErrorCategory{};
}  // namespace

inline std::error_code make_error_code(MicmRosenbrockErrc e)
{
  return { static_cast<int>(e), rosenbrockErrorCategory };
}

namespace micm
{

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
      case SolverState::InfDetected: return "InfDetected";
      default: return "Unknown";
    }
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::RosenbrockSolver()
      : processes_(),
        parameters_(RosenbrockSolverParameters::three_stage_rosenbrock_parameters()),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
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
            [](const std::vector<Process>& processes,
               const std::map<std::string, std::size_t>& variable_map) -> ProcessSetPolicy {
              return ProcessSetPolicy{ processes, variable_map };
            })
  {
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::RosenbrockSolver(
      const System& system,
      const std::vector<Process>& processes,
      const RosenbrockSolverParameters& parameters,
      const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
      const std::function<ProcessSetPolicy(const std::vector<Process>&, const std::map<std::string, std::size_t>&)>
          create_process_set)
      : processes_(processes),
        parameters_(parameters),
        state_parameters_(),
        process_set_(),
        linear_solver_()
  {
    MICM_PROFILE_FUNCTION();

    std::map<std::string, std::size_t> variable_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;

    std::size_t index = 0;
    auto used_species = ProcessSetPolicy::SpeciesUsed(processes);
    auto available_species = system.UniqueNames();
    std::sort(available_species.begin(), available_species.end());
    std::set<std::string> unused_species;
    std::set_difference(
        available_species.begin(),
        available_species.end(),
        used_species.begin(),
        used_species.end(),
        std::inserter(unused_species, unused_species.begin()));
    if (unused_species.size() > 0 && !parameters_.ignore_unused_species_)
    {
      std::string err_msg = "Unused species in chemical system:";
      for (auto& species : unused_species)
        err_msg += " '" + species + "'";
      err_msg += ".";
      throw std::system_error(make_error_code(MicmRosenbrockErrc::UnusedSpecies), err_msg);
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

    // if the tolerances aren't already set, initialize them and then set based off of information in the system
    if (parameters_.absolute_tolerance_.size() != variable_map.size())
    {
      parameters_.absolute_tolerance_ = std::vector<double>(variable_map.size(), 1e-3);
      for (auto& species : system.gas_phase_.species_)
      {
        if (species.HasProperty("absolute tolerance"))
        {
          parameters_.absolute_tolerance_[variable_map[species.name_]] = species.GetProperty<double>("absolute tolerance");
        }
      }
      for (auto& phase : system.phases_)
      {
        for (auto& species : phase.second.species_)
        {
          if (species.HasProperty("absolute tolerance"))
          {
            parameters_.absolute_tolerance_[variable_map[species.name_]] = species.GetProperty<double>("absolute tolerance");
          }
        }
      }
    }

    // setup the state_parameters
    std::vector<std::string> param_labels{};
    for (const auto& process : processes_)
      if (process.rate_constant_)
        for (auto& label : process.rate_constant_->CustomParameters())
          param_labels.push_back(label);

    process_set_ = std::move(create_process_set(processes, variable_map));

    auto jacobian = build_jacobian<SparseMatrixPolicy>(
        process_set_.NonZeroJacobianElements(), parameters_.number_of_grid_cells_, system.StateSize());

    std::vector<std::size_t> jacobian_diagonal_elements;
    jacobian_diagonal_elements.reserve(jacobian[0].size());
    for (std::size_t i = 0; i < jacobian[0].size(); ++i)
      jacobian_diagonal_elements.push_back(jacobian.VectorIndex(0, i, i));

    state_parameters_ = { .number_of_grid_cells_ = parameters_.number_of_grid_cells_,
                          .number_of_species_ = variable_map.size(),
                          .number_of_rate_constants_ = processes_.size(),
                          .variable_names_ = system.UniqueNames(state_reordering),
                          .custom_rate_parameter_labels_ = param_labels,
                          .jacobian_diagonal_elements_ = jacobian_diagonal_elements,
                          .nonzero_jacobian_elements_ = process_set_.NonZeroJacobianElements() };

    process_set_.SetJacobianFlatIds(jacobian);
    linear_solver_ = std::move(create_linear_solver(jacobian, 1.0e-30));
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy>
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::GetState() const
  {
    return State<MatrixPolicy, SparseMatrixPolicy>{ state_parameters_ };
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::SolverResult
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::Solve(
      double time_step,
      State<MatrixPolicy, SparseMatrixPolicy>& state) noexcept
  {
    MICM_PROFILE_FUNCTION();

    typename RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::SolverResult result{};
    result.state_ = SolverState::Running;
    // reset the upper, lower matrix. Repeated calls without zeroing these matrices can lead to lack of convergence
    auto& lower = state.lower_matrix_.AsVector();
    auto& upper = state.upper_matrix_.AsVector();
    for (auto& l : lower)
      l = 0;
    for (auto& u : upper)
      u = 0;
    MatrixPolicy<double> Y(state.variables_);
    std::size_t num_rows = Y.NumRows();
    std::size_t num_cols = Y.NumColumns();
    MatrixPolicy<double> Ynew(num_rows, num_cols, 0.0);
    MatrixPolicy<double> initial_forcing(num_rows, num_cols, 0.0);
    MatrixPolicy<double> forcing(num_rows, num_cols, 0.0);
    MatrixPolicy<double> temp(num_rows, num_cols, 0.0);
    std::vector<MatrixPolicy<double>> K{};
    const double h_max = parameters_.h_max_ == 0.0 ? time_step : std::min(time_step, parameters_.h_max_);
    const double h_start =
        parameters_.h_start_ == 0.0 ? std::max(parameters_.h_min_, delta_min_) : std::min(h_max, parameters_.h_start_);

    SolverStats stats;

    UpdateState(state);

    K.reserve(parameters_.stages_);
    for (std::size_t i = 0; i < parameters_.stages_; ++i)
      K.emplace_back(num_rows, num_cols, 0.0);

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
      CalculateForcing(state.rate_constants_, Y, initial_forcing);
      stats.function_calls += 1;

      // compute the negative jacobian at the beginning of the current time
      CalculateNegativeJacobian(state.rate_constants_, Y, state.jacobian_);
      stats.jacobian_updates += 1;

      bool accepted = false;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        bool is_singular{ false };
        // Form and factor the rosenbrock ode jacobian
        LinearFactor(H, parameters_.gamma_[0], is_singular, Y, stats, state);
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
                Ynew.Axpy(parameters_.a_[stage_combinations + j], K[j]);
              }
              CalculateForcing(state.rate_constants_, Ynew, forcing);
              stats.function_calls += 1;
            }
          }
          K[stage].AsVector().assign(forcing.AsVector().begin(), forcing.AsVector().end());
          for (uint64_t j = 0; j < stage; ++j)
          {
            K[stage].Axpy(parameters_.c_[stage_combinations + j] / H, K[j]);
          }
          temp.AsVector().assign(K[stage].AsVector().begin(), K[stage].AsVector().end());
          linear_solver_.template Solve<MatrixPolicy>(temp, K[stage], state.lower_matrix_, state.upper_matrix_);
          stats.solves += 1;
        }

        // Compute the new solution
        Ynew.AsVector().assign(Y.AsVector().begin(), Y.AsVector().end());
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Ynew.Axpy(parameters_.m_[stage], K[stage]);

        // Compute the error estimation
        MatrixPolicy<double> Yerror(num_rows, num_cols, 0);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Yerror.Axpy(parameters_.e_[stage], K[stage]);

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
        else if (std::isinf(error) == 1)
        {
          Y.AsVector().assign(Ynew.AsVector().begin(), Ynew.AsVector().end());
          result.state_ = SolverState::InfDetected;
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

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing)
  {
    MICM_PROFILE_FUNCTION();
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
    process_set_.template AddForcingTerms<MatrixPolicy>(rate_constants, number_densities, forcing);
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const requires(!VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    MICM_PROFILE_FUNCTION();

    for (std::size_t i_block = 0; i_block < jacobian.Size(); ++i_block)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_block * jacobian.FlatBlockSize());
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::AlphaMinusJacobian(
      SparseMatrixPolicy<double>& jacobian,
      const double& alpha) const requires(VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    MICM_PROFILE_FUNCTION();

    const std::size_t n_cells = jacobian.GroupVectorSize();
    for (std::size_t i_group = 0; i_group < jacobian.NumberOfGroups(jacobian.Size()); ++i_group)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_group * jacobian.GroupSize(jacobian.FlatBlockSize()));
      for (const auto& i_elem : state_parameters_.jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void
  RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::CalculateNegativeJacobian(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      SparseMatrixPolicy<double>& jacobian)
  {
    MICM_PROFILE_FUNCTION();

    std::fill(jacobian.AsVector().begin(), jacobian.AsVector().end(), 0.0);
    process_set_.template SubtractJacobianTerms<MatrixPolicy, SparseMatrixPolicy>(
        rate_constants, number_densities, jacobian);
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::UpdateState(
      State<MatrixPolicy, SparseMatrixPolicy>& state)
  {
    Process::UpdateState(processes_, state);
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline void RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::LinearFactor(
      double& H,
      const double gamma,
      bool& singular,
      const MatrixPolicy<double>& number_densities,
      SolverStats& stats,
      State<MatrixPolicy, SparseMatrixPolicy>& state)
  {
    MICM_PROFILE_FUNCTION();

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
      }
      else
      {
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

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline double RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::NormalizedError(
      const MatrixPolicy<double>& Y,
      const MatrixPolicy<double>& Ynew,
      const MatrixPolicy<double>& errors) const requires(!VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    MICM_PROFILE_FUNCTION();

    auto& _y = Y.AsVector();
    auto& _ynew = Ynew.AsVector();
    auto& _errors = errors.AsVector();
    std::size_t N = Y.AsVector().size();
    std::size_t n_species = state_parameters_.number_of_species_;

    double ymax = 0;
    double errors_over_scale = 0;
    double error = 0;

    for (std::size_t i = 0; i < N; ++i)
    {
      ymax = std::max(std::abs(_y[i]), std::abs(_ynew[i]));
      errors_over_scale =
          _errors[i] / (parameters_.absolute_tolerance_[i % n_species] + parameters_.relative_tolerance_ * ymax);
      error += errors_over_scale * errors_over_scale;
    }

    double error_min = 1.0e-10;

    return std::max(std::sqrt(error / N), error_min);
  }

  template<
      template<class>
      class MatrixPolicy,
      template<class>
      class SparseMatrixPolicy,
      class LinearSolverPolicy,
      class ProcessSetPolicy>
  inline double RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>::NormalizedError(
      const MatrixPolicy<double>& Y,
      const MatrixPolicy<double>& Ynew,
      const MatrixPolicy<double>& errors) const requires(VectorizableSparse<SparseMatrixPolicy<double>>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    MICM_PROFILE_FUNCTION();

    auto y_iter = Y.AsVector().begin();
    auto ynew_iter = Ynew.AsVector().begin();
    auto errors_iter = errors.AsVector().begin();
    std::size_t N = Y.NumRows() * Y.NumColumns();
    const std::size_t L = Y.GroupVectorSize();
    std::size_t n_species = state_parameters_.number_of_species_;

    std::size_t whole_blocks = std::floor(Y.NumRows() / Y.GroupVectorSize()) * Y.GroupSize();

    double errors_over_scale = 0;
    double error = 0;

    // compute the error over the blocks which fit exactly into the L parameter
    for (std::size_t i = 0; i < whole_blocks; ++i)
    {
      errors_over_scale =
          *errors_iter / (parameters_.absolute_tolerance_[(i / L) % n_species] +
                          parameters_.relative_tolerance_ * std::max(std::abs(*y_iter), std::abs(*ynew_iter)));
      error += errors_over_scale * errors_over_scale;
      ++y_iter;
      ++ynew_iter;
      ++errors_iter;
    }

    // compute the error over the remaining elements that are in the next group but didn't fill a full group
    std::size_t remaining_rows = Y.NumRows() % Y.GroupVectorSize();

    if (remaining_rows > 0)
    {
      for (std::size_t y = 0; y < Y.NumColumns(); ++y)
      {
        for (std::size_t x = 0; x < remaining_rows; ++x)
        {
          std::size_t idx = y * L + x;
          errors_over_scale = errors_iter[idx] /
                              (parameters_.absolute_tolerance_[y] +
                               parameters_.relative_tolerance_ * std::max(std::abs(y_iter[idx]), std::abs(ynew_iter[idx])));
          error += errors_over_scale * errors_over_scale;
        }
      }
    }

    double error_min = 1.0e-10;
    return std::max(std::sqrt(error / N), error_min);
  }

}  // namespace micm
