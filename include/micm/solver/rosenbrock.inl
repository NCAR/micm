// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/reducers.hpp>
#include <micm/util/types.hpp>

namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class StatePolicy>
  inline SolverResult AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::Solve(
      Real time_step,
      StatePolicy& state,
      const RosenbrockSolverParameters& parameters) const noexcept
  {
    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);

    SolverResult result{};
    result.state_ = SolverState::Running;
    DenseMatrixPolicy& Y = state.variables_;  // Y will hold the new solution at the end of the solve
    auto derived_class_temporary_variables =
        static_cast<RosenbrockTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    DenseMatrixPolicy& Ynew = derived_class_temporary_variables->Ynew_;
    DenseMatrixPolicy& initial_forcing = derived_class_temporary_variables->initial_forcing_;
    std::vector<DenseMatrixPolicy>& K = derived_class_temporary_variables->K_;
    DenseMatrixPolicy& Yerror = derived_class_temporary_variables->Yerror_;
    auto& current_c_over_h = derived_class_temporary_variables->current_c_over_h_;
    auto& error = derived_class_temporary_variables->error_;
    auto& diagonal = state.views_.upper_left_identity_diagonal_;

    const Real h_min = parameters.h_min_ == 0.0 ? DEFAULT_H_MIN * time_step : parameters.h_min_;
    const Real h_max = parameters.h_max_ == 0.0 ? time_step : std::min(time_step, parameters.h_max_);
    const Real h_start = parameters.h_start_ == 0.0 ? DEFAULT_H_START * time_step : std::min(h_max, parameters.h_start_);
    Real H = std::min(std::max(h_min, std::abs(h_start)), std::abs(h_max));

    const bool has_constraints = constraints_.Size() > 0;

    // Declared here so they remain in scope for the solver loop below (captured by reference by mass_coupling).
    // std::function gives mass_coupling a concrete, nameable type; the closure type returned by
    // DenseMatrixPolicy::Function() is anonymous and cannot be named directly.
    std::function<void(DenseMatrixPolicy&, DenseMatrixPolicy&)> mass_coupling;

    // Initialize algebraic constraint variables and pre-build the mass-coupling function.
    // All K matrices have the same shape, so K[0] is used to capture column counts at creation time.
    if (has_constraints)
    {
      mass_coupling = DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& k_stage_view,
              const typename DenseMatrixPolicy::ConstViewType& k_j_view) {
            for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
            {
              if (diagonal[i_var] != 0.0)
              {
                k_stage_view.ForEachRow(
                    [&current_c_over_h](Real& ks, const Real& kj) { ks += current_c_over_h * kj; },
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

    Real present_time = 0.0;

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
      Real last_alpha = 0.0;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        // Compute alpha for AlphaMinusJacobian function
        Real alpha = 1.0 / (H * parameters.gamma_[0]);
        if constexpr (!LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          // The Jacobian retains the alpha shift applied on earlier attempts of
          // this step, so shift by the difference only. last_alpha must hold the
          // cumulative shift now present in the matrix, not the per-attempt delta.
          const double cumulative_alpha = alpha;
          alpha -= last_alpha;
          last_alpha = cumulative_alpha;
        }

        // Form and factor the rosenbrock ode jacobian
        LinearFactor(alpha, result.stats_, state);

        // Compute the stages
        for (Index stage = 0; stage < parameters.stages_; ++stage)
        {
          const Index stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
          if (stage == 0)
          {
            K[stage].Copy(initial_forcing);
          }
          else
          {
            if (parameters.new_function_evaluation_[stage])
            {
              Ynew.Copy(Y);
              for (Index j = 0; j < stage; ++j)
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
          for (Index j = 0; j < stage; ++j)
          {
            Real c_over_h = parameters.c_[stage_combinations + j] / H;
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
              current_c_over_h.CopyToDevice();
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
        for (Index stage = 0; stage < parameters.stages_; ++stage)
        {
          Ynew.Axpy(parameters.m_[stage], K[stage]);
        }

        Yerror.Fill(0);
        for (Index stage = 0; stage < parameters.stages_; ++stage)
        {
          Yerror.Axpy(parameters.e_[stage], K[stage]);
        }

        // For DAE systems, replace the near-zero algebraic error estimates with step changes.
        // The embedded error formula produces Yerror ≈ 0 for algebraic rows (M_ii = 0 zeroes
        // the inter-stage coupling terms), making the solver insensitive to algebraic tolerances.
        // Setting Yerror[a] = Ynew[a] - Y[a] allows the solver to reject steps where algebraic
        // variables change more than their tolerance permits, preventing overshoot.
        //
        // IMPORTANT: This uses the step change as-is without order scaling. For very stiff systems
        // where the embedded error estimate K[3] ≈ 0 for all variables (including differential),
        // this provides the ONLY non-trivial error signal for algebraic variables. The step change
        // is O(H) while the true error is O(H^(p+1)), so this is conservative — the solver may
        // take more steps than strictly necessary. However, it correctly prevents overshoot by
        // rejecting steps where algebraic variables change more than their tolerance permits.
        if (has_constraints)
        {
          constraints_.SetAlgebraicErrors(Yerror, Y, Ynew);
        }

        // Compute the normalized error
        static_cast<const Derived*>(this)->NormalizedError(Y, Ynew, Yerror, state, error);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
        Real fac = std::min(
            parameters.factor_max_,
            std::max(
                parameters.factor_min_,
                parameters.safety_factor_ / std::pow(error, 1 / parameters.estimator_of_local_order_)));
        Real Hnew = H * fac;

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

    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class SparseMatrixPolicy, class StatePolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::AlphaMinusJacobian(
      StatePolicy& state,
      const Real& alpha) const
  {
    // Form [alpha * M - J] by scaling diagonal updates with the mass matrix diagonal.
    // ODE rows have M[i][i]=1 and get +alpha; algebraic rows have M[i][i]=0 and get no alpha shift.
    auto& views = state.views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
          Index i_diag = 0;
          for (const auto& i_elem : views.jacobian_diagonal_elements_)
          {
            const Real scaled_alpha = alpha * views.upper_left_identity_diagonal_[i_diag++];
            jacobian_view.ForEachBlock(
                [scaled_alpha](Real& diag) { diag += scaled_alpha; }, jacobian_view.GetBlockView(i_elem));
          }
        },
        state.jacobian_)(state.jacobian_);
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class StatePolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::LinearFactor(
      const Real alpha,
      SolverStats& stats,
      StatePolicy& state) const
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
  template<class DenseMatrixPolicy, class StatePolicy>
  inline void AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::NormalizedError(
      const DenseMatrixPolicy& Y,
      const DenseMatrixPolicy& Ynew,
      const DenseMatrixPolicy& errors,
      const StatePolicy& state,
      typename DenseMatrixPolicy::template ScalarType<Real>& error) const
  {
    using SumType = typename DenseMatrixPolicy::template SumType<Real>;

    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7
    const auto& atol = std::as_const(state.absolute_tolerance_).GetView();
    const auto& rtol = state.relative_tolerance_;
    const Index n_vars = Y.NumColumns();
    const Index n_cells = Y.NumRows();

    error = 0;
    error.CopyToDevice();

    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ConstViewType& y_view,
            const typename DenseMatrixPolicy::ConstViewType& ynew_view,
            const typename DenseMatrixPolicy::ConstViewType& errors_view) {
          for (Index i_var = 0; i_var < n_vars; ++i_var)
          {
            // skip padding rows so their possibly non-zero values
            // do not end up in the normalized error.
            y_view.ReduceStrict(
                SumType{ error },
                [&](const Real& y, const Real& ynew, const Real& var_error, Real& acc)
                {
                  Real ymax = (std::abs(y) > std::abs(ynew) ? std::abs(y) : std::abs(ynew));
                  Real errors_over_scale = var_error / (atol[i_var % n_vars] + rtol * ymax);
                  acc += errors_over_scale * errors_over_scale;
                },
                y_view.GetConstColumnView(i_var),
                ynew_view.GetConstColumnView(i_var),
                errors_view.GetConstColumnView(i_var));
          }
        },
        Y,
        Ynew,
        errors)(Y, Ynew, errors);
    error.CopyToHost();
    constexpr Real error_min = 1.0e-10;
    const Index N = std::max<Index>(1, Y.NumRows() * Y.NumColumns());

    error = std::max(std::sqrt(error / N), error_min);
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
  template<class StatePolicy>
  inline SolverState
  AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy, Derived>::InitializeConstraints(
      StatePolicy& state,
      const RosenbrockSolverParameters& parameters,
      SolverStats& stats) const noexcept
  {
    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);
    using LOrType = typename DenseMatrixPolicy::LOrType;
    using MaxType = typename DenseMatrixPolicy::template MaxType<Real>;

    auto& Y = state.variables_;
    auto& diagonal = state.views_.upper_left_identity_diagonal_;
    auto derived_class_temporary_variables =
        static_cast<RosenbrockTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    // Reuse initial_forcing_ as the residual/delta workspace, K_[0] as the rollback buffer, Ynew_ as
    // the line-search candidate, and Yerror_ as the candidate's simplified Newton correction: no
    // stage vector is read, and neither Ynew_ nor Yerror_ is read, until the first stage of the
    // integration loop below.
    auto& delta = derived_class_temporary_variables->initial_forcing_;
    auto& original_variables = derived_class_temporary_variables->K_[0];
    auto& candidate = derived_class_temporary_variables->Ynew_;
    auto& candidate_correction = derived_class_temporary_variables->Yerror_;
    auto& max_residual = derived_class_temporary_variables->max_residual_;
    auto& max_correction = derived_class_temporary_variables->max_correction_;
    auto& nan_detected = derived_class_temporary_variables->nan_detected_;
    auto& inf_detected = derived_class_temporary_variables->inf_detected_;
    auto& step = derived_class_temporary_variables->constraint_init_step_;
    const auto& atol = std::as_const(state.absolute_tolerance_).GetView();
    const auto& rtol = state.relative_tolerance_;
    const Index n_vars = Y.NumColumns();

    // Pre-build reusable Function objects outside the iteration loop
    // ||.||_inf over the algebraic rows only, with non-finite values flagged
    auto check_algebraic_values = DenseMatrixPolicy::Function(
        MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& delta_view) {
          for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
          {
            if (diagonal[i_var] == 0.0)
            {
              auto col_view = delta_view.GetConstColumnView(i_var);
              delta_view.ReduceStrict(
                  LOrType{ nan_detected }, [](const Real& val, Bool& acc) { acc = acc || std::isnan(val); }, col_view);

              delta_view.ReduceStrict(
                  LOrType{ inf_detected }, [](const Real& val, Bool& acc) { acc = acc || std::isinf(val); }, col_view);

              // exclude padded cells incase they are non-zero
              delta_view.ReduceStrict(
                  MaxType{ max_residual },
                  [](const Real& val, Real& acc)
                  {
                    Real abs_val = std::abs(val);
                    if (!std::isnan(abs_val) && !std::isinf(abs_val))
                    {
                      acc = (acc > abs_val ? acc : abs_val);
                    }
                  },
                  col_view);
            }
          }
        },
        delta);

    // Measure of the Newton correction in the same weighted state space that defines
    // integration accuracy:
    //   scale_a = atol_a + rtol * max(|z_a|, |z_a + delta_a|)
    //   q = max_a |delta_a| / scale_a
    // q is dimensionless, and scaling a complete constraint row leaves it unchanged.
    auto check_weighted_correction = DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ConstViewType& y_view,
            const typename DenseMatrixPolicy::ConstViewType& delta_view) {
          for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
          {
            if (diagonal[i_var] == 0.0)
            {
              // exclude padded cells incase they are non-zero
              y_view.ReduceStrict(
                  MaxType{ max_correction },
                  [&](const Real& z, const Real& correction, Real& acc)
                  {
                    Real corrected = z + correction;
                    Real z_max = (std::abs(z) > std::abs(corrected) ? std::abs(z) : std::abs(corrected));
                    Real scale = atol[i_var % n_vars] + rtol * z_max;
                    Real q = scale > 0.0 ? std::abs(correction) / scale
                                         : (correction == 0.0 ? 0.0 : std::numeric_limits<Real>::infinity());
                    acc = (acc > q ? acc : q);
                  },
                  y_view.GetConstColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    auto apply_update = DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& y_view,
            const typename DenseMatrixPolicy::ConstViewType& delta_view) {
          for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
          {
            if (diagonal[i_var] == 0.0)
            {
              auto d_col_view = delta_view.GetConstColumnView(i_var);
              y_view.ReduceStrict(
                  LOrType{ nan_detected }, [](const Real& d_val, Bool& acc) { acc = acc || std::isnan(d_val); }, d_col_view);
              y_view.ReduceStrict(
                  LOrType{ inf_detected }, [](const Real& d_val, Bool& acc) { acc = acc || std::isinf(d_val); }, d_col_view);
              y_view.ForEachRow(
                  [&](Real& y_val, const Real& d_val)
                  {
                    if (!std::isnan(d_val) && !std::isinf(d_val))
                    {
                      y_val += d_val;
                    }
                  },
                  y_view.GetColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    // Line-search candidate: candidate += step * delta on the algebraic rows only. Identical to
    // apply_update above apart from the step factor, so a damped candidate and an undamped update
    // treat padded rows and non-finite corrections the same way. The step length has to be a
    // device-resident scalar rather than a captured Real: MICM_LAMBDA captures by value, so a Real
    // named here would be frozen at the value it held when this Function object was built. Same
    // mechanism the stage loop uses for current_c_over_h_.
    auto set_candidate = DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& candidate_view,
            const typename DenseMatrixPolicy::ConstViewType& delta_view) {
          for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
          {
            if (diagonal[i_var] == 0.0)
            {
              candidate_view.ForEachRow(
                  [&](Real& c_val, const Real& d_val)
                  {
                    if (!std::isnan(d_val) && !std::isinf(d_val))
                    {
                      c_val += step * d_val;
                    }
                  },
                  candidate_view.GetColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    // Projection is transactional: a failure must not hand the caller a partially updated
    // algebraic state, which is neither their input nor a solution. The snapshot is taken
    // lazily, so a projection that converges without applying an update costs no copy.
    bool variables_modified = false;
    auto restore_and_return = [&](SolverState status)
    {
      if (variables_modified)
      {
        Y.Copy(original_variables);
      }
      return status;
    };

    // One pass beyond the update limit measures the correction that remains after the final
    // permitted update, so a state that reaches the manifold on that update is not reported
    // as a failure.
    for (Index update = 0; update <= parameters.constraint_init_max_iterations_; ++update)
    {
      // 1. Evaluate constraint residuals: G(y)
      delta.Fill(0);
      constraints_.AddForcingTerms(Y, state.custom_rate_parameters_, delta);

      // 2. Reject non-finite residuals. A residual has whatever units and scale its constraint
      //    equation was written in, so it cannot decide convergence; only an exactly zero
      //    residual is already consistent and is accepted without a factorization.
      max_residual = 0;
      nan_detected = false;
      inf_detected = false;
      max_residual.CopyToDevice();
      nan_detected.CopyToDevice();
      inf_detected.CopyToDevice();
      check_algebraic_values(delta);
      max_residual.CopyToHost();
      nan_detected.CopyToHost();
      inf_detected.CopyToHost();

      if (nan_detected)
      {
        return restore_and_return(SolverState::NaNDetected);
      }
      if (inf_detected)
      {
        return restore_and_return(SolverState::InfDetected);
      }

      stats.constraint_init_iterations_ += 1;

      if (max_residual == 0.0)
      {
        return SolverState::Converged;
      }

      // 3. Compute constraint Jacobian: -dG/dy
      state.jacobian_.Fill(0);
      constraints_.SubtractJacobianTerms(Y, state.custom_rate_parameters_, state.jacobian_);
      stats.jacobian_updates_ += 1;

      // 4. Form system matrix with alpha=1.0:
      //    ODE rows (M=1): identity on diagonal (constraint Jacobian has no ODE-row entries)
      //    Algebraic rows (M=0): -dG/dy from SubtractJacobianTerms (no alpha contribution)
      static_cast<const Derived*>(this)->template AlphaMinusJacobian<SparseMatrixPolicy>(state, 1.0);

      // 5. Factor the system matrix
      if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
      {
        linear_solver_.Factor(state.jacobian_);
      }
      else
      {
        linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
      }
      stats.decompositions_ += 1;

      // 6. Solve: -J * delta = G  =>  delta = -J^{-1} * G
      if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
      {
        linear_solver_.Solve(delta, state.jacobian_);
      }
      else
      {
        linear_solver_.Solve(delta, state.lower_matrix_, state.upper_matrix_);
      }
      stats.solves_ += 1;

      // 7. Reject a non-finite Newton correction before it is measured or applied
      max_residual = 0;
      nan_detected = false;
      inf_detected = false;
      max_residual.CopyToDevice();
      nan_detected.CopyToDevice();
      inf_detected.CopyToDevice();
      check_algebraic_values(delta);
      max_residual.CopyToHost();
      nan_detected.CopyToHost();
      inf_detected.CopyToHost();

      if (nan_detected)
      {
        return restore_and_return(SolverState::NaNDetected);
      }
      if (inf_detected)
      {
        return restore_and_return(SolverState::InfDetected);
      }

      // 8. Converged once the remaining correction is a small fraction of the state-error
      //    scale set by the configured absolute and relative tolerances
      max_correction = 0;
      max_correction.CopyToDevice();
      check_weighted_correction(Y, delta);
      max_correction.CopyToHost();

      if (max_correction <= parameters.constraint_init_tolerance_)
      {
        return SolverState::Converged;
      }

      if (update == parameters.constraint_init_max_iterations_)
      {
        break;
      }

      // 9. Apply update only to algebraic variables, snapshotting the caller's state first
      if (!variables_modified)
      {
        original_variables.Copy(Y);
        variables_modified = true;
      }

      if (parameters.constraint_init_max_backtracks_ == 0)
      {
        // Line search disabled. Every statement in this block is unchanged from the undamped
        // implementation, in the same order on the same data, so this path is bit-for-bit identical.
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        apply_update(Y, delta);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected)
        {
          return restore_and_return(SolverState::NaNDetected);
        }
        if (inf_detected)
        {
          return restore_and_return(SolverState::InfDetected);
        }
        continue;
      }

      // 10. Damped update: a backtracking line search on the step length. The merit is the weighted
      //     norm of the simplified Newton correction at the candidate -- the correction the Jacobian
      //     already factored at step 5 produces for the candidate's residual. Scaling a complete
      //     constraint row scales G and dG/dy together, so it leaves both corrections and therefore
      //     every decision below unchanged, exactly as it leaves the step-8 test unchanged. Reusing
      //     the factorization keeps a trial to one constraint forcing evaluation and one triangular
      //     solve.
      //
      //     Step 8 has already reduced ||delta(Y)||_w into max_correction, so the sufficient-decrease
      //     reference value costs no extra kernel and no extra device round trip. It is strictly
      //     greater than constraint_init_tolerance_ here, which keeps the threshold below away from
      //     denormals in any precision.
      const Real current_correction_norm = max_correction;

      Real step_length = 1.0;
      bool update_accepted = false;
      for (Index backtrack = 0; backtrack <= parameters.constraint_init_max_backtracks_; ++backtrack)
      {
        // 10a. Candidate iterate y + step * delta, algebraic rows only. Step 7 established that
        //      delta is finite there, so a non-finite candidate means the sum overflowed. Reject it
        //      here rather than downstream: the equilibrium residual clamps negative concentrations
        //      to zero and would report a finite residual for an infinite candidate.
        step = step_length;
        step.CopyToDevice();
        candidate.Copy(Y);
        set_candidate(candidate, delta);

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // 10b. Candidate residual G(y + step * delta)
        candidate_correction.Fill(0);
        constraints_.AddForcingTerms(candidate, state.custom_rate_parameters_, candidate_correction);

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate_correction);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // An exactly zero candidate residual is already consistent, by the same rule step 2 applies
        // at Y. Accepting it here also keeps a singular factorization from turning an exactly zero
        // right-hand side into 0/0 in the solve below.
        if (max_residual == 0.0)
        {
          update_accepted = true;
          break;
        }

        // 10c. Simplified Newton correction, reusing the factorization from step 5. Solve reads the
        //      factors through a const view and never writes them, so one factorization serves every
        //      trial.
        if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          linear_solver_.Solve(candidate_correction, state.jacobian_);
        }
        else
        {
          linear_solver_.Solve(candidate_correction, state.lower_matrix_, state.upper_matrix_);
        }
        stats.solves_ += 1;

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate_correction);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // 10d. Accept when the candidate's remaining correction is converged, or when it decreases
        //      the correction at Y by the required fraction.
        max_correction = 0;
        max_correction.CopyToDevice();
        check_weighted_correction(candidate, candidate_correction);
        max_correction.CopyToHost();

        const Real candidate_correction_norm = max_correction;
        const Real required_norm =
            (Real{ 1.0 } - parameters.constraint_init_sufficient_decrease_ * step_length) * current_correction_norm;
        if (candidate_correction_norm <= parameters.constraint_init_tolerance_ ||
            candidate_correction_norm < required_norm)
        {
          update_accepted = true;
          break;
        }

        step_length *= parameters.constraint_init_backtrack_factor_;
      }

      if (!update_accepted)
      {
        return restore_and_return(SolverState::ConstraintInitializationFailed);
      }

      // Take the accepted candidate whole, so the state carried into the next update is exactly the
      // one whose merit was measured. Convergence is deliberately not reported from inside the line
      // search: the next pass re-evaluates the residual at the new state and decides through the
      // step-2 and step-8 tests, which keeps the meaning of constraint_init_iterations_ unchanged.
      Y.Copy(candidate);
    }

    // Did not converge within the permitted number of Newton updates
    return restore_and_return(SolverState::ConstraintInitializationFailed);
  }

}  // namespace micm
