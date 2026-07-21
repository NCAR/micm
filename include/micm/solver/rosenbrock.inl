// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  inline SolverResult AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::Solve(
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
    if (parameters.h_persist_ && state.solver_step_size_suggestion_ > 0.0)
    {
      H = std::min(std::max(h_min, state.solver_step_size_suggestion_), std::abs(h_max));
    }
    // Controller suggestion carried into the next Solve() call when h_persist_ is
    // set; tracked separately so end-of-interval truncation cannot poison it.
    double h_persist_suggestion = H;

    const bool has_constraints = constraints_.Size() > 0;

    // Declared here so they remain in scope for the solver loop below (captured by reference by mass_coupling).
    // std::function gives mass_coupling a concrete, nameable type; the closure type returned by
    // DenseMatrixPolicy::Function() is anonymous and cannot be named directly.
    double current_c_over_h = 0.0;
    const auto& diagonal = state.upper_left_identity_diagonal_;
    std::function<void(DenseMatrixPolicy&, DenseMatrixPolicy&)> mass_coupling;

    // Initialize algebraic constraint variables and pre-build the mass-coupling function.
    // All K matrices have the same shape, so K[0] is used to capture column counts at creation time.
    if (has_constraints)
    {
      mass_coupling = DenseMatrixPolicy::Function(
          [&current_c_over_h, &diagonal](auto&& k_stage_view, auto&& k_j_view)
          {
            for (std::size_t i_var = 0; i_var < diagonal.size(); ++i_var)
            {
              if (diagonal[i_var] != 0.0)
              {
                k_stage_view.ForEachRow(
                    [&current_c_over_h](double& ks, const double& kj) { ks += current_c_over_h * kj; },
                    k_stage_view.GetColumnView(i_var),
                    k_j_view.GetConstColumnView(i_var));
              }
            }
          },
          K[0],
          K[0]);

      auto init_state = InitializeConstraints(state, parameters, result.stats_);
      if (init_state != SolverState::Converged)
      {
        result.state_ = init_state;
        return result;
      }
    }

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
      const double h_before_end_clamp = H;
      H = std::min(H, std::abs(time_step - present_time));

      // compute the initial forcing at the beginning of the current time
      initial_forcing.Fill(0);
      rates_.AddForcingTerms(state, Y, initial_forcing);

      if (has_constraints)
      {
        constraints_.AddForcingTerms(Y, state.custom_rate_parameters_, initial_forcing);
      }

      result.stats_.function_calls_ += 1;

      // compute the negative jacobian at the beginning of the current time
      state.jacobian_.Fill(0);
      rates_.SubtractJacobianTerms(state, Y, state.jacobian_);

      if (has_constraints)
      {
        constraints_.SubtractJacobianTerms(Y, state.custom_rate_parameters_, state.jacobian_);
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
          const std::size_t stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
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
              rates_.AddForcingTerms(state, Ynew, K[stage]);
              if (has_constraints)
              {
                constraints_.AddForcingTerms(Ynew, state.custom_rate_parameters_, K[stage]);
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
            if (!has_constraints)
            {
              K[stage].Axpy(c_over_h, K[j]);
            }
            else
            {
              // DAE path: scale by mass matrix diagonal element-wise.
              // For ODE variables (diagonal = 1), accumulate c/H * K[j].
              // For algebraic variables (diagonal = 0), the coupling is zero.
              current_c_over_h = c_over_h;
              mass_coupling(K[stage], K[j]);
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

        Ynew.Copy(Y);
        for (uint64_t stage = 0; stage < parameters.stages_; ++stage)
        {
          Ynew.Axpy(parameters.m_[stage], K[stage]);
        }

        // Embedded local truncation error estimate, formed over ALL rows
        // (differential and algebraic). This is the method's true O(H^(p+1))
        // error estimate. For index-1 DAE systems the algebraic rows correctly
        // evaluate to ~0: an algebraic variable is slaved to the differential
        // variables through its constraint, so it carries negligible independent
        // local error — its accuracy follows from the differential variables
        // (which are error-controlled here) plus constraint enforcement at each
        // stage. We therefore use this estimate as-is for all variables.
        //
        // NOTE: a previous implementation overwrote the algebraic rows with the
        // raw step change Ynew[a]-Y[a]. That is O(H), not a truncation error, so
        // dividing it by a tight algebraic tolerance throttled the step size to
        // hold dY_algebraic ~ atol per step — inflating step counts (up to ~100x
        // on stiff problems) for no accuracy gain, while providing no feasibility
        // protection that the constraint solve does not already provide. See
        // docs/superpowers/notes/2026-05-28-dae-algebraic-error-investigation.md.
        Yerror.Fill(0);
        for (uint64_t stage = 0; stage < parameters.stages_; ++stage)
        {
          Yerror.Axpy(parameters.e_[stage], K[stage]);
        }

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
        if (std::isinf(error) == 1)
        {
          result.state_ = SolverState::InfDetected;
          break;
        }
        if ((error < 1) || (H < h_min))
        {
          result.stats_.accepted_ += 1;
          present_time = present_time + H;
          Y.Swap(Ynew);
          // A step shortened only to land on the interval end says nothing about
          // stability at the controller's natural size, so the pre-clamp value is
          // the better continuation suggestion; a rejection within this segment
          // invalidates it.
          const bool truncated_only = (H < h_before_end_clamp) && !reject_last_h;
          Hnew = std::max(h_min, std::min(Hnew, h_max));
          if (reject_last_h)
          {
            // No step size increase after a rejected step
            Hnew = std::min(Hnew, H);
          }
          reject_last_h = false;
          reject_more_h = false;
          H = Hnew;
          h_persist_suggestion = truncated_only ? std::max(Hnew, std::min(h_before_end_clamp, h_max)) : Hnew;
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
            rates_.SubtractJacobianTerms(state, Y, state.jacobian_);
            // Subtract constraint Jacobian terms (for DAE systems)
            if (has_constraints)
            {
              constraints_.SubtractJacobianTerms(Y, state.custom_rate_parameters_, state.jacobian_);
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

    if (parameters.h_persist_ && result.state_ == SolverState::Converged)
    {
      state.solver_step_size_suggestion_ = h_persist_suggestion;
    }

    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::AlphaMinusJacobian(
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

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class SparseMatrixPolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::AlphaMinusJacobian(
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
        {
          jacobian_vector[i_elem + i_cell] += alpha * diagonal_scale;
        }
      }
    }
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::LinearFactor(
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

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class DenseMatrixPolicy>
  inline double AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::NormalizedError(
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
    const std::size_t N = std::max<std::size_t>(1, Y.NumRows() * Y.NumColumns());

    return std::max(std::sqrt(error / N), error_min);
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class DenseMatrixPolicy>
  inline double AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::NormalizedError(
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
    const std::size_t N = std::max<std::size_t>(1, Y.NumRows() * Y.NumColumns());

    return std::max(std::sqrt(error / N), error_min);
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  inline SolverState
  AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::InitializeConstraints(
      auto& state,
      const RosenbrockSolverParameters& parameters,
      SolverStats& stats) const noexcept
  {
    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);

    auto& Y = state.variables_;
    auto derived_class_temporary_variables =
        static_cast<RosenbrockTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    auto& correction = derived_class_temporary_variables->initial_forcing_;
    auto& candidate = derived_class_temporary_variables->Ynew_;
    auto& candidate_correction = derived_class_temporary_variables->Yerror_;
    auto& original_state = derived_class_temporary_variables->K_[0];

    const auto& diagonal = state.upper_left_identity_diagonal_;
    const auto& atol = state.absolute_tolerance_;
    const double rtol = state.relative_tolerance_;
    const std::size_t number_of_tolerances = atol.size();

    // Constraint initialization is transactional: a failed projection must not leave a partially
    // updated algebraic state for the caller to mistake for a valid solution.
    original_state.Copy(Y);

    auto finite_algebraic_values = [&diagonal](const DenseMatrixPolicy& values)
    {
      SolverState result = SolverState::Converged;
      for (std::size_t i_cell = 0; i_cell < values.NumRows(); ++i_cell)
      {
        for (std::size_t i_var = 0; i_var < diagonal.size(); ++i_var)
        {
          if (diagonal[i_var] != 0.0)
          {
            continue;
          }
          if (std::isnan(values[i_cell][i_var]))
          {
            return SolverState::NaNDetected;
          }
          if (std::isinf(values[i_cell][i_var]))
          {
            result = SolverState::InfDetected;
          }
        }
      }
      return result;
    };

    auto maximum_algebraic_magnitude = [&diagonal](const DenseMatrixPolicy& values)
    {
      double maximum = 0.0;
      for (std::size_t i_cell = 0; i_cell < values.NumRows(); ++i_cell)
      {
        for (std::size_t i_var = 0; i_var < diagonal.size(); ++i_var)
        {
          if (diagonal[i_var] == 0.0)
          {
            maximum = std::max(maximum, std::abs(values[i_cell][i_var]));
          }
        }
      }
      return maximum;
    };

    auto weighted_correction_norm =
        [&diagonal, &atol, rtol, number_of_tolerances](const DenseMatrixPolicy& values, const DenseMatrixPolicy& updates)
    {
      double maximum = 0.0;
      for (std::size_t i_cell = 0; i_cell < values.NumRows(); ++i_cell)
      {
        for (std::size_t i_var = 0; i_var < diagonal.size(); ++i_var)
        {
          if (diagonal[i_var] != 0.0)
          {
            continue;
          }
          const double value = values[i_cell][i_var];
          const double update = updates[i_cell][i_var];
          const double scale =
              atol[i_var % number_of_tolerances] + rtol * std::max(std::abs(value), std::abs(value + update));
          const double normalized_update =
              scale > 0.0 ? std::abs(update) / scale : (update == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
          maximum = std::max(maximum, normalized_update);
        }
      }
      return maximum;
    };

    auto set_candidate =
        [&diagonal](
            DenseMatrixPolicy& destination, const DenseMatrixPolicy& values, const DenseMatrixPolicy& updates, double step)
    {
      destination.Copy(values);
      for (std::size_t i_cell = 0; i_cell < values.NumRows(); ++i_cell)
      {
        for (std::size_t i_var = 0; i_var < diagonal.size(); ++i_var)
        {
          if (diagonal[i_var] == 0.0)
          {
            destination[i_cell][i_var] += step * updates[i_cell][i_var];
          }
        }
      }
    };

    auto restore_and_return = [&Y, &original_state](SolverState result)
    {
      Y.Copy(original_state);
      return result;
    };

    for (std::size_t iter = 0; iter < parameters.constraint_init_max_iterations_; ++iter)
    {
      // Evaluate G(y). Exact zeros are accepted without a factorization; all non-zero
      // residuals are judged through their state-space Newton correction below.
      correction.Fill(0);
      constraints_.AddForcingTerms(Y, state.custom_rate_parameters_, correction);
      stats.constraint_init_iterations_ += 1;
      const auto residual_state = finite_algebraic_values(correction);
      if (residual_state != SolverState::Converged)
      {
        return restore_and_return(residual_state);
      }
      if (maximum_algebraic_magnitude(correction) == 0.0)
      {
        return SolverState::Converged;
      }

      // Compute and factor -dG/dy with identity rows for differential variables.
      state.jacobian_.Fill(0);
      constraints_.SubtractJacobianTerms(Y, state.custom_rate_parameters_, state.jacobian_);
      stats.jacobian_updates_ += 1;
      static_cast<const Derived*>(this)->template AlphaMinusJacobian<SparseMatrixPolicy>(state, 1.0);
      if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
      {
        linear_solver_.Factor(state.jacobian_);
      }
      else
      {
        linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
      }
      stats.decompositions_ += 1;

      // Solve -J delta = G, so delta is the Newton correction in state units.
      if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
      {
        linear_solver_.Solve(correction, state.jacobian_);
      }
      else
      {
        linear_solver_.Solve(correction, state.lower_matrix_, state.upper_matrix_);
      }
      stats.solves_ += 1;

      const auto correction_state = finite_algebraic_values(correction);
      if (correction_state != SolverState::Converged)
      {
        return restore_and_return(correction_state);
      }

      const double current_correction_norm = weighted_correction_norm(Y, correction);
      if (!std::isfinite(current_correction_norm))
      {
        return restore_and_return(SolverState::InfDetected);
      }

      // Globalize Newton with an affine-covariant merit function: the norm of the
      // simplified Newton correction computed with the already-factored current
      // Jacobian. Scaling a complete constraint row changes neither correction.
      double step = 1.0;
      bool update_accepted = false;
      double accepted_correction_norm = std::numeric_limits<double>::infinity();
      for (std::size_t backtrack = 0; backtrack <= parameters.constraint_init_max_backtracks_; ++backtrack)
      {
        set_candidate(candidate, Y, correction, step);
        if (finite_algebraic_values(candidate) != SolverState::Converged)
        {
          step *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        candidate_correction.Fill(0);
        constraints_.AddForcingTerms(candidate, state.custom_rate_parameters_, candidate_correction);
        if (finite_algebraic_values(candidate_correction) != SolverState::Converged)
        {
          step *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        if (maximum_algebraic_magnitude(candidate_correction) == 0.0)
        {
          accepted_correction_norm = 0.0;
          update_accepted = true;
          break;
        }

        if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          linear_solver_.Solve(candidate_correction, state.jacobian_);
        }
        else
        {
          linear_solver_.Solve(candidate_correction, state.lower_matrix_, state.upper_matrix_);
        }
        stats.solves_ += 1;

        if (finite_algebraic_values(candidate_correction) != SolverState::Converged)
        {
          step *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        accepted_correction_norm = weighted_correction_norm(candidate, candidate_correction);
        const double required_norm =
            (1.0 - parameters.constraint_init_sufficient_decrease_ * step) * current_correction_norm;
        if (accepted_correction_norm <= parameters.constraint_init_tolerance_ || accepted_correction_norm < required_norm)
        {
          update_accepted = true;
          break;
        }

        step *= parameters.constraint_init_backtrack_factor_;
      }

      if (!update_accepted)
      {
        return restore_and_return(SolverState::ConstraintInitializationFailed);
      }

      Y.Copy(candidate);
      if (accepted_correction_norm <= parameters.constraint_init_tolerance_)
      {
        return SolverState::Converged;
      }
    }

    return restore_and_return(SolverState::ConstraintInitializationFailed);
  }

}  // namespace micm
