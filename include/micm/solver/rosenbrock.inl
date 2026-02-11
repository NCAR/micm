// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
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

    const bool has_constraints = constraints_.Size() > 0;

    double present_time = 0.0;

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_step + parameters.round_off_) <= 0 && (result.state_ == SolverState::Running))
    {
      if (result.stats_.number_of_steps_ > parameters.max_number_of_steps_)
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

      // Add constraint residuals to forcing (for DAE systems)
      if (has_constraints)
      {
        constraints_.AddForcingTerms(Y, initial_forcing);
      }
      result.stats_.function_calls_ += 1;

      // compute the negative jacobian at the beginning of the current time
      state.jacobian_.Fill(0);
      rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);

      // Add constraint Jacobian terms (for DAE systems)
      if (has_constraints)
      {
        constraints_.SubtractJacobianTerms(Y, state.jacobian_);
      }
      result.stats_.jacobian_updates_ += 1;

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
        LinearFactor(alpha, result.stats_, state);

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
              // Copy state variables from Y to Ynew for the new function evaluation
              Ynew.Fill(0.0);
              for (std::size_t i_cell = 0; i_cell < Y.NumRows(); ++i_cell)
              {
                for (std::size_t i_var = 0; i_var < state.state_size_; ++i_var)
                {
                  Ynew[i_cell][i_var] = Y[i_cell][i_var];
                }
              }
              for (uint64_t j = 0; j < stage; ++j)
              {
                Ynew.Axpy(parameters.a_[stage_combinations + j], K[j]);
              }
              K[stage].Fill(0);
              rates_.AddForcingTerms(state.rate_constants_, Ynew, K[stage]);
              // Add constraint residuals for DAE systems
              if (has_constraints)
              {
                constraints_.AddForcingTerms(Ynew, K[stage]);
              }
              result.stats_.function_calls_ += 1;
            }
          }
          if (stage + 1 < parameters.stages_ && !parameters.new_function_evaluation_[stage + 1])
          {
            K[stage + 1].Copy(K[stage]);
          }
          for (uint64_t j = 0; j < stage; ++j)
          {
            const double c_over_h = parameters.c_[stage_combinations + j] / H;
            for (std::size_t i_cell = 0; i_cell < K[stage].NumRows(); ++i_cell)
            {
              for (std::size_t i_var = 0; i_var < K[stage].NumColumns(); ++i_var)
              {
                K[stage][i_cell][i_var] +=
                    c_over_h * state.upper_left_identity_diagonal_[i_var] * K[j][i_cell][i_var];
              }
            }
          }
          if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
          {
            linear_solver_.Solve(K[stage], state.jacobian_);
          }
          else
          {
            linear_solver_.Solve(K[stage], state.lower_matrix_, state.upper_matrix_);
          }
          result.stats_.solves_ += 1;
        }

        // Compute the new solution: Copy Y, then add increments from K stages
        for (std::size_t i_cell = 0; i_cell < Y.NumRows(); ++i_cell)
        {
          for (std::size_t i_var = 0; i_var < state.state_size_; ++i_var)
          {
            Ynew[i_cell][i_var] = Y[i_cell][i_var];
          }
        }
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

        result.stats_.number_of_steps_ += 1;

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
          result.stats_.accepted_ += 1;
          present_time = present_time + H;
          // Copy solution from Ynew back to Y
          for (std::size_t i_cell = 0; i_cell < Y.NumRows(); ++i_cell)
          {
            for (std::size_t i_var = 0; i_var < state.state_size_; ++i_var)
            {
              Y[i_cell][i_var] = Ynew[i_cell][i_var];
            }
          }
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
          if (result.stats_.accepted_ >= 1)
          {
            result.stats_.rejected_ += 1;
          }
          // Re-generate the Jacobian matrix for the inline LU algorithm
          if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
          {
            state.jacobian_.Fill(0);
            rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);
            // Subtract constraint Jacobian terms (for DAE systems)
            if (has_constraints)
            {
              constraints_.SubtractJacobianTerms(Y, state.jacobian_);
            }
            result.stats_.jacobian_updates_ += 1;
          }
        }
      }
    }

    if (result.state_ == SolverState::Running)
    {
      result.state_ = SolverState::Converged;
    }

    result.stats_.final_time_ = present_time;

    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, Derived>::AlphaMinusJacobian(
      auto& state,
      const double& alpha) const
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  {
    // Form [alpha * M - J] by scaling diagonal updates with the mass matrix diagonal.
    // ODE rows have M[i][i]=1 and get +alpha; algebraic rows have M[i][i]=0 and get no alpha shift.
    for (std::size_t i_block = 0; i_block < state.jacobian_.NumberOfBlocks(); ++i_block)
    {
      auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_block * state.jacobian_.FlatBlockSize());
      std::size_t i_diag = 0;
      for (const auto& i_elem : state.jacobian_diagonal_elements_)
      {
        jacobian_vector[i_elem] += alpha * state.upper_left_identity_diagonal_[i_diag++];
      }
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
    // Form [alpha * M - J] by scaling diagonal updates with the mass matrix diagonal.
    for (std::size_t i_group = 0; i_group < state.jacobian_.NumberOfGroups(state.jacobian_.NumberOfBlocks()); ++i_group)
    {
      auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_group * state.jacobian_.GroupSize());
      std::size_t i_diag = 0;
      for (const auto& i_elem : state.jacobian_diagonal_elements_)
      {
        const double diagonal_scale = state.upper_left_identity_diagonal_[i_diag++];
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          jacobian_vector[i_elem + i_cell] += alpha * diagonal_scale;
      }
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

    const auto& atol = state.absolute_tolerance_;
    const auto& rtol = state.relative_tolerance_;
    const std::size_t N = Y.NumRows() * Y.NumColumns();
    const std::size_t n_vars = atol.size();

    double ymax = 0;
    double errors_over_scale = 0;
    double error = 0;

    for (std::size_t i_cell = 0; i_cell < Y.NumRows(); ++i_cell)
    {
      for (std::size_t i_var = 0; i_var < Y.NumColumns(); ++i_var)
      {
        ymax = std::max(std::abs(Y[i_cell][i_var]), std::abs(Ynew[i_cell][i_var]));
        errors_over_scale = errors[i_cell][i_var] / (atol[i_var % n_vars] + rtol * ymax);
        error += errors_over_scale * errors_over_scale;
      }
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

    const auto& atol = state.absolute_tolerance_;
    const auto& rtol = state.relative_tolerance_;
    const std::size_t N = Y.NumRows() * Y.NumColumns();
    const std::size_t n_vars = atol.size();

    double ymax = 0;
    double errors_over_scale = 0;
    double error = 0;

    // Use row/column indexing so error estimation is independent of matrix storage layout.
    for (std::size_t i_cell = 0; i_cell < Y.NumRows(); ++i_cell)
    {
      for (std::size_t i_var = 0; i_var < Y.NumColumns(); ++i_var)
      {
        ymax = std::max(std::abs(Y[i_cell][i_var]), std::abs(Ynew[i_cell][i_var]));
        errors_over_scale = errors[i_cell][i_var] / (atol[i_var % n_vars] + rtol * ymax);
        error += errors_over_scale * errors_over_scale;
      }
    }

    double error_min = 1.0e-10;
    return std::max(std::sqrt(error / N), error_min);
  }

}  // namespace micm
