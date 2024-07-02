// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  inline SolverResult AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::Solve(
      double time_step,
      auto& state) noexcept
  {
    MICM_PROFILE_FUNCTION();
    using MatrixPolicy = decltype(state.variables_);

    SolverResult result{};
    result.state_ = SolverState::Running;
    MatrixPolicy Y = state.variables_;
    std::size_t num_rows = Y.NumRows();
    std::size_t num_cols = Y.NumColumns();
    MatrixPolicy Ynew(num_rows, num_cols, 0.0);
    MatrixPolicy initial_forcing(num_rows, num_cols, 0.0);
    MatrixPolicy forcing(num_rows, num_cols, 0.0);
    std::vector<MatrixPolicy> K{};
    const double h_max = parameters_.h_max_ == 0.0 ? time_step : std::min(time_step, parameters_.h_max_);
    const double h_start =
        parameters_.h_start_ == 0.0 ? std::max(parameters_.h_min_, DELTA_MIN) : std::min(h_max, parameters_.h_start_);

    SolverStats stats;

    K.reserve(parameters_.stages_);
    for (std::size_t i = 0; i < parameters_.stages_; ++i)
      K.emplace_back(num_rows, num_cols, 0.0);

    double present_time = 0.0;
    double H = std::min(std::max(std::abs(parameters_.h_min_), std::abs(h_start)), std::abs(h_max));

    if (std::abs(H) <= 10 * parameters_.round_off_)
      H = DELTA_MIN;

    // TODO: the logic above this point should be moved to the constructor and should return an error
    //       if the parameters are invalid (e.g., h_min > h_max or h_start > h_max)

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_step + parameters_.round_off_) <= 0 && (result.state_ == SolverState::Running))
    {
      if (stats.number_of_steps_ > parameters_.max_number_of_steps_)
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
      initial_forcing.Fill(0.0);
      rates_.AddForcingTerms(state.rate_constants_, Y, initial_forcing);
      stats.function_calls_ += 1;

      // compute the negative jacobian at the beginning of the current time
      state.jacobian_.Fill(0.0);
      rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);
      stats.jacobian_updates_ += 1;

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
              Ynew.Copy(Y);
              for (uint64_t j = 0; j < stage; ++j)
              {
                Ynew.Axpy(parameters_.a_[stage_combinations + j], K[j]);
              }
              forcing.Fill(0.0);
              rates_.AddForcingTerms(state.rate_constants_, Ynew, forcing);
              stats.function_calls_ += 1;
            }
          }
          K[stage].Copy(forcing);
          for (uint64_t j = 0; j < stage; ++j)
          {
            K[stage].Axpy(parameters_.c_[stage_combinations + j] / H, K[j]);
          }
          linear_solver_.Solve(K[stage], state.lower_matrix_, state.upper_matrix_);
          stats.solves_ += 1;
        }

        // Compute the new solution
        Ynew.Copy(Y);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Ynew.Axpy(parameters_.m_[stage], K[stage]);

        // Compute the error estimation
        MatrixPolicy Yerror(num_rows, num_cols, 0);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
          Yerror.Axpy(parameters_.e_[stage], K[stage]);

        // Compute the normalized error
        auto error = static_cast<Derived*>(this)->NormalizedError(Y, Ynew, Yerror);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
        double fac = std::min(
            parameters_.factor_max_,
            std::max(
                parameters_.factor_min_,
                parameters_.safety_factor_ / std::pow(error, 1 / parameters_.estimator_of_local_order_)));
        double Hnew = H * fac;

        stats.number_of_steps_ += 1;

        // Check the error magnitude and adjust step size
        if (std::isnan(error))
        {
          Y.Copy(Ynew);
          result.state_ = SolverState::NaNDetected;
          break;
        }
        else if (std::isinf(error) == 1)
        {
          Y.Copy(Ynew);
          result.state_ = SolverState::InfDetected;
          break;
        }
        else if ((error < 1) || (H < parameters_.h_min_))
        {
          stats.accepted_ += 1;
          present_time = present_time + H;
          Y.Copy(Ynew);
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
          if (stats.accepted_ >= 1)
          {
            stats.rejected_ += 1;
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
    state.variables_ = Y;

    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::AlphaMinusJacobian(
      SparseMatrixPolicy& jacobian,
      const double& alpha) const requires(!VectorizableSparse<SparseMatrixPolicy>)
  {
    MICM_PROFILE_FUNCTION();

    for (std::size_t i_block = 0; i_block < jacobian.NumberOfBlocks(); ++i_block)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_block * jacobian.FlatBlockSize());
      for (const auto& i_elem : jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::AlphaMinusJacobian(
      SparseMatrixPolicy& jacobian,
      const double& alpha) const requires(VectorizableSparse<SparseMatrixPolicy>)
  {
    MICM_PROFILE_FUNCTION();

    const std::size_t n_cells = jacobian.GroupVectorSize();
    for (std::size_t i_group = 0; i_group < jacobian.NumberOfGroups(jacobian.NumberOfBlocks()); ++i_group)
    {
      auto jacobian_vector = std::next(jacobian.AsVector().begin(), i_group * jacobian.GroupSize(jacobian.FlatBlockSize()));
      for (const auto& i_elem : jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::LinearFactor(
      double& H,
      const double gamma,
      bool& singular,
      const auto& number_densities,
      SolverStats& stats,
      auto& state)
  {
    MICM_PROFILE_FUNCTION();

    auto jacobian = state.jacobian_;
    uint64_t n_consecutive = 0;
    singular = false;
    while (true)
    {
      double alpha = 1 / (H * gamma);
      static_cast<Derived*>(this)->AlphaMinusJacobian(jacobian, alpha);
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
      stats.decompositions_ += 1;
      if (!singular)
        break;
      stats.singular_ += 1;
      if (++n_consecutive > 5)
        break;
      H /= 2;
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class DenseMatrixPolicy>
  inline double AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::NormalizedError(
      const DenseMatrixPolicy& Y,
      const DenseMatrixPolicy& Ynew,
      const DenseMatrixPolicy& errors) const requires(!VectorizableDense<DenseMatrixPolicy>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    MICM_PROFILE_FUNCTION();

    auto& _y = Y.AsVector();
    auto& _ynew = Ynew.AsVector();
    auto& _errors = errors.AsVector();
    std::size_t N = Y.AsVector().size();
    std::size_t n_vars = parameters_.absolute_tolerance_.size();

    double ymax = 0;
    double errors_over_scale = 0;
    double error = 0;

    for (std::size_t i = 0; i < N; ++i)
    {
      ymax = std::max(std::abs(_y[i]), std::abs(_ynew[i]));
      errors_over_scale =
          _errors[i] / (parameters_.absolute_tolerance_[i % n_vars] + parameters_.relative_tolerance_ * ymax);
      error += errors_over_scale * errors_over_scale;
    }

    double error_min = 1.0e-10;

    return std::max(std::sqrt(error / N), error_min);
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class DenseMatrixPolicy>
  inline double AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::NormalizedError(
      const DenseMatrixPolicy& Y,
      const DenseMatrixPolicy& Ynew,
      const DenseMatrixPolicy& errors) const requires(VectorizableDense<DenseMatrixPolicy>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    MICM_PROFILE_FUNCTION();

    auto y_iter = Y.AsVector().begin();
    auto ynew_iter = Ynew.AsVector().begin();
    auto errors_iter = errors.AsVector().begin();
    std::size_t N = Y.NumRows() * Y.NumColumns();
    const std::size_t L = Y.GroupVectorSize();
    std::size_t n_vars = parameters_.absolute_tolerance_.size();

    std::size_t whole_blocks = std::floor(Y.NumRows() / Y.GroupVectorSize()) * Y.GroupSize();

    double errors_over_scale = 0;
    double error = 0;

    // compute the error over the blocks which fit exactly into the L parameter
    for (std::size_t i = 0; i < whole_blocks; ++i)
    {
      errors_over_scale =
          *errors_iter / (parameters_.absolute_tolerance_[(i / L) % n_vars] +
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
