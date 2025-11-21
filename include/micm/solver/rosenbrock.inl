// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  inline SolverResult AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::Solve(
      double time_step,
      auto& state,
      const RosenbrockSolverParameters& parameters) const noexcept
  {
    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);

    SolverResult result{};
    result.state_ = SolverState::Running;
    auto& Y = state.variables_;  // Y will hold the new solution at the end of the solve
    auto derived_class_temporary_variables =
        static_cast<RosenbrockTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    auto& Ynew = derived_class_temporary_variables->Ynew_;
    auto& initial_forcing = derived_class_temporary_variables->initial_forcing_;
    auto& K = derived_class_temporary_variables->K_;
    auto& Yerror = derived_class_temporary_variables->Yerror_;
    const double h_min = parameters.h_min_ == 0.0 ? DEFAULT_H_MIN * time_step : parameters.h_min_;
    const double h_max = parameters.h_max_ == 0.0 ? time_step : std::min(time_step, parameters.h_max_);
    const double h_start = parameters.h_start_ == 0.0 ? DEFAULT_H_START * time_step : std::min(h_max, parameters.h_start_);
    double H = std::min(std::max(h_min, std::abs(h_start)), std::abs(h_max));

    SolverStats stats;

    double present_time = 0.0;

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_step + parameters.round_off_) <= 0 && (result.state_ == SolverState::Running))
    {
      if (stats.number_of_steps_ > parameters.max_number_of_steps_)
      {
        result.state_ = SolverState::ConvergenceExceededMaxSteps;
        break;
      }

      if (((present_time + 0.1 * H) == present_time) || (H <= parameters.round_off_))
      {
        result.state_ = SolverState::StepSizeTooSmall;
        break;
      }

      //  Limit H if necessary to avoid going beyond the specified chemistry time step
      H = std::min(H, std::abs(time_step - present_time));

      // compute the initial forcing at the beginning of the current time
      initial_forcing.Fill(0);
      rates_.AddForcingTerms(state.rate_constants_, Y, initial_forcing);
      stats.function_calls_ += 1;

      // compute the negative jacobian at the beginning of the current time
      state.jacobian_.Fill(0);
      rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);
      stats.jacobian_updates_ += 1;

      bool accepted = false;
      double last_alpha = 0.0;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        // Compute alpha for AlphaMinusJacobian function
        double alpha = 1.0 / (H * parameters.gamma_[0]);
        if constexpr (!LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          // Compute alpha accounting for the last alpha value
          // This is necessary to avoid the need to re-factor the jacobian for non-inline LU algorithms
          alpha -= last_alpha;
          last_alpha = alpha;
        }

        // Form and factor the rosenbrock ode jacobian
        LinearFactor(alpha, stats, state);

        // Compute the stages
        for (uint64_t stage = 0; stage < parameters.stages_; ++stage)
        {
          double stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
          if (stage == 0)
          {
            K[stage].Copy(initial_forcing);
          }
          else
          {
            if (parameters.new_function_evaluation_[stage])
            {
              Ynew.Copy(Y);
              for (uint64_t j = 0; j < stage; ++j)
              {
                Ynew.Axpy(parameters.a_[stage_combinations + j], K[j]);
              }
              K[stage].Fill(0);
              rates_.AddForcingTerms(state.rate_constants_, Ynew, K[stage]);
              stats.function_calls_ += 1;
            }
          }
          if (stage + 1 < parameters.stages_ && !parameters.new_function_evaluation_[stage + 1])
          {
            K[stage + 1].Copy(K[stage]);
          }
          for (uint64_t j = 0; j < stage; ++j)
          {
            K[stage].Axpy(parameters.c_[stage_combinations + j] / H, K[j]);
          }
          if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
          {
            linear_solver_.Solve(K[stage], state.jacobian_);
          }
          else
          {
            linear_solver_.Solve(K[stage], state.lower_matrix_, state.upper_matrix_);
          }
          stats.solves_ += 1;
        }

        // Compute the new solution
        Ynew.Copy(Y);
        for (uint64_t stage = 0; stage < parameters.stages_; ++stage)
          Ynew.Axpy(parameters.m_[stage], K[stage]);

        Yerror.Fill(0);
        for (uint64_t stage = 0; stage < parameters.stages_; ++stage)
          Yerror.Axpy(parameters.e_[stage], K[stage]);

        // Compute the normalized error
        auto error = static_cast<const Derived*>(this)->NormalizedError(Y, Ynew, Yerror, state);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
        double fac = std::min(
            parameters.factor_max_,
            std::max(
                parameters.factor_min_,
                parameters.safety_factor_ / std::pow(error, 1 / parameters.estimator_of_local_order_)));
        double Hnew = H * fac;

        stats.number_of_steps_ += 1;

        // Check the error magnitude and adjust step size
        if (std::isnan(error))
        {
          result.state_ = SolverState::NaNDetected;
          break;
        }
        else if (std::isinf(error) == 1)
        {
          result.state_ = SolverState::InfDetected;
          break;
        }
        else if ((error < 1) || (H < h_min))
        {
          stats.accepted_ += 1;
          present_time = present_time + H;
          Y.Swap(Ynew);
          Hnew = std::max(h_min, std::min(Hnew, h_max));
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
            Hnew = H * parameters.rejection_factor_decrease_;
          }
          reject_more_h = reject_last_h;
          reject_last_h = true;
          H = Hnew;
          if (stats.accepted_ >= 1)
          {
            stats.rejected_ += 1;
          }
          // Re-generate the Jacobian matrix for the inline LU algorithm
          if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
          {
            state.jacobian_.Fill(0);
            rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);
            stats.jacobian_updates_ += 1;
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

    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::AlphaMinusJacobian(
      auto& state,
      const double& alpha) const
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  {
    for (std::size_t i_block = 0; i_block < state.jacobian_.NumberOfBlocks(); ++i_block)
    {
      auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_block * state.jacobian_.FlatBlockSize());
      for (const auto& i_elem : state.jacobian_diagonal_elements_)
        jacobian_vector[i_elem] += alpha;
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::AlphaMinusJacobian(
      auto& state,
      const double& alpha) const
    requires(VectorizableSparse<SparseMatrixPolicy>)
  {
    constexpr std::size_t n_cells = SparseMatrixPolicy::GroupVectorSize();
    for (std::size_t i_group = 0; i_group < state.jacobian_.NumberOfGroups(state.jacobian_.NumberOfBlocks()); ++i_group)
    {
      auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_group * state.jacobian_.GroupSize());
      for (const auto& i_elem : state.jacobian_diagonal_elements_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha;
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::LinearFactor(
      const double alpha,
      SolverStats& stats,
      auto& state) const
  {
    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);

    static_cast<const Derived*>(this)->template AlphaMinusJacobian<SparseMatrixPolicy>(state, alpha);

    if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
    {
      linear_solver_.Factor(state.jacobian_);
    }
    else
    {
      linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
    }
    stats.decompositions_ += 1;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class DenseMatrixPolicy>
  inline double AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::NormalizedError(
      const DenseMatrixPolicy& Y,
      const DenseMatrixPolicy& Ynew,
      const DenseMatrixPolicy& errors,
      auto& state) const
    requires(!VectorizableDense<DenseMatrixPolicy>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    auto& _y = Y.AsVector();
    auto& _ynew = Ynew.AsVector();
    auto& _errors = errors.AsVector();
    const auto& atol = state.absolute_tolerance_;
    const auto& rtol = state.relative_tolerance_;
    const std::size_t N = Y.AsVector().size();
    const std::size_t n_vars = atol.size();

    double ymax = 0;
    double errors_over_scale = 0;
    double error = 0;

    for (std::size_t i = 0; i < N; ++i)
    {
      ymax = std::max(std::abs(_y[i]), std::abs(_ynew[i]));
      errors_over_scale = _errors[i] / (atol[i % n_vars] + rtol * ymax);
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
      const DenseMatrixPolicy& errors,
      auto& state) const
    requires(VectorizableDense<DenseMatrixPolicy>)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7

    auto y_iter = Y.AsVector().begin();
    auto ynew_iter = Ynew.AsVector().begin();
    auto errors_iter = errors.AsVector().begin();
    const auto& atol = state.absolute_tolerance_;
    auto rtol = state.relative_tolerance_;
    const std::size_t N = Y.NumRows() * Y.NumColumns();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    const std::size_t whole_blocks = std::floor(Y.NumRows() / Y.GroupVectorSize()) * Y.GroupSize();
    const std::size_t n_vars = atol.size();

    double errors_over_scale = 0;
    double error = 0;

    // compute the error over the blocks which fit exactly into the L parameter
    for (std::size_t i = 0; i < whole_blocks; ++i)
    {
      errors_over_scale = *errors_iter / (atol[(i / L) % n_vars] + rtol * std::max(std::abs(*y_iter), std::abs(*ynew_iter)));
      error += errors_over_scale * errors_over_scale;
      ++y_iter;
      ++ynew_iter;
      ++errors_iter;
    }

    // compute the error over the remaining elements that are in the next group but didn't fill a full group
    const std::size_t remaining_rows = Y.NumRows() % Y.GroupVectorSize();

    if (remaining_rows > 0)
    {
      for (std::size_t y = 0; y < Y.NumColumns(); ++y)
      {
        for (std::size_t x = 0; x < remaining_rows; ++x)
        {
          const std::size_t idx = y * L + x;
          errors_over_scale =
              errors_iter[idx] / (atol[y] + rtol * std::max(std::abs(y_iter[idx]), std::abs(ynew_iter[idx])));
          error += errors_over_scale * errors_over_scale;
        }
      }
    }

    double error_min = 1.0e-10;
    return std::max(std::sqrt(error / N), error_min);
  }

}  // namespace micm
