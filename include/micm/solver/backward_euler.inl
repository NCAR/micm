// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/reducers.hpp>
#include <micm/util/types.hpp>

namespace micm
{
  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
  inline SolverResult AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>::Solve(
      Real time_step,
      auto& state,
      const BackwardEulerSolverParameters& parameters) const
  {
    // A fully implicit euler implementation is given by the following equation:
    // y_{n+1} = y_n + H * f(t_{n+1}, y_{n+1})
    // This is a root finding problem because you need to know y_{n+1} to compute f(t_{n+1}, y_{n+1})
    // you need to solve the equation y_{n+1} - y_n - H f(t_{n+1}, y_{n+1}) = 0
    // A series of time step reductions are used after failed solves to try to find a solution
    // These reductions are controlled by the time_step_reductions parameter in the solver parameters
    // if the last attempt to reduce the timestep fails,
    // accept the current H but do not update the Yn vector

    using DenseMatrixPolicy = decltype(state.variables_);
    using SparseMatrixPolicy = decltype(state.jacobian_);

    SolverResult result;

    Index max_iter = parameters.max_number_of_steps_;
    const auto time_step_reductions = parameters.time_step_reductions_;

    Real H = parameters.h_start_ == 0.0 ? time_step : parameters.h_start_;
    Real present_time = 0.0;
    Index n_successful_integrations = 0;
    Index n_convergence_failures = 0;

    auto derived_class_temporary_variables =
        static_cast<BackwardEulerTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    auto& Yn = derived_class_temporary_variables->Yn_;
    auto& Yn1 = state.variables_;  // Yn1 will hold the new solution at the end of the solve
    auto& forcing = derived_class_temporary_variables->forcing_;

    // Ensure Yn starts with the same values as the state variables
    Yn.Copy(Yn1);

    while (present_time < time_step)
    {
      result.state_ = SolverState::Running;
      bool converged = false;
      Index iterations = 0;

      do
      {
        result.stats_.number_of_steps_++;
        // the first time Yn1 is equal to Yn
        // after the first iteration Yn1 is updated to the new solution
        // so we can use Yn1 to calculate the forcing and jacobian
        // calculate forcing
        forcing.Fill(0.0);
        rates_.AddForcingTerms(state, Yn1, forcing);
        result.stats_.function_calls_++;

        // calculate the negative jacobian
        state.jacobian_.Fill(0.0);
        rates_.SubtractJacobianTerms(state, Yn1, state.jacobian_);
        result.stats_.jacobian_updates_++;

        // add the inverse of the time step from the diagonal
        state.jacobian_.AddToDiagonal(1 / H);

        // We want to solve this equation for a zero
        // (y_{n+1} - y_n) / H = f(t_{n+1}, y_{n+1})

        // try to find the root by factoring and solving the linear system
        if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          linear_solver_.Factor(state.jacobian_);
        }
        else
        {
          linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
        }
        result.stats_.decompositions_++;

        // forcing_blk in camchem
        // residual = forcing - (Yn1 - Yn) / H
        // since forcing is only used once, we can reuse it to store the residual
        forcing.ForEach([&](Real& f, const Real& yn1, const Real& yn) { f -= (yn1 - yn) / H; }, Yn1, Yn);

        // the result of the linear solver will be stored in forcing
        // this represents the change in the solution
        if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          linear_solver_.Solve(forcing, state.jacobian_);
        }
        else
        {
          linear_solver_.Solve(forcing, state.lower_matrix_, state.upper_matrix_);
        }
        result.stats_.solves_++;

        // solution_blk in camchem
        // Yn1 = Yn1 + residual;
        // always make sure the solution is positive regardless of which iteration we are on
        Yn1.ForEach([&](Real& yn1, const Real& f) { yn1 = std::max<Real>(0.0, yn1 + f); }, forcing);

        // if this is the first iteration, we don't need to check for convergence
        if (iterations++ == 0)
        {
          continue;
        }

        // check for convergence
        converged = IsConverged(parameters, forcing, Yn1, state.absolute_tolerance_, state.relative_tolerance_);
      } while (!converged && iterations < max_iter);

      if (!converged)
      {
        result.stats_.rejected_++;
        n_successful_integrations = 0;

        if (n_convergence_failures >= time_step_reductions.size())
        {
          present_time += H;
          // Distinguish a genuine blow-up from a merely unconverged step, so callers are not handed
          // non-finite concentrations under a "success" status. Only reached once the solver has
          // already exhausted its step-size reductions, so this scan costs nothing in the happy path.
          result.state_ = SolverState::AcceptingUnconvergedIntegration;
          for (const auto& y : Yn1.AsVector())
          {
            if (std::isnan(y))
            {
              result.state_ = SolverState::NaNDetected;
              break;
            }
            if (std::isinf(y))
            {
              result.state_ = SolverState::InfDetected;
            }
          }
          break;
        }

        // if we fail, we need to reset the solution to the last known good solution
        Yn1.Copy(Yn);
        H *= time_step_reductions[n_convergence_failures++];
      }
      else
      {
        result.stats_.accepted_++;
        result.state_ = SolverState::Converged;
        present_time += H;
        Yn.Copy(Yn1);

        // when we accept two solutions in a row, we can increase the time step
        n_successful_integrations++;
        if (n_successful_integrations >= 2)
        {
          n_successful_integrations = 0;
          H *= 2.0;
        }
      }
      // Don't let H go past the time step
      H = std::min(H, time_step - present_time);
    }

    result.stats_.final_time_ = present_time;
    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
  template<class DenseMatrixPolicy>
  inline bool AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>::IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1,
      const std::vector<Real>& absolute_tolerance,
      Real relative_tolerance)
  {
    const Index n_vars = absolute_tolerance.size();
    const Real small = parameters.small_;
    const Real rel_tol = relative_tolerance;
    bool retval = true;
    DenseMatrixPolicy::Function(
        [&](const auto&& residual_view, const auto&& Yn1_view)
        {
          for (Index i_var = 0; i_var < n_vars; ++i_var)
          {
            const Real var_abs_tol = absolute_tolerance[i_var];
            residual_view.ReduceStrict(
                LAnd{ retval },
                [=](const Real& residual, const Real& Yn1, bool& acc)
                {
                  // A non-finite residual is never converged. Without this check an infinite residual escapes
                  // the test below, because the relative bound rel_tol * |Yn1| is itself infinite and
                  // inf > inf is false -- a blown-up solve would then be reported as SolverState::Converged.
                  if (!std::isfinite(residual))
                  {
                    acc = false;
                    return;
                  }
                  if (std::abs(residual) > small && std::abs(residual) > var_abs_tol &&
                      std::abs(residual) > rel_tol * std::abs(Yn1))
                  {
                    acc = false;
                  }
                },
                residual_view.GetConstColumnView(i_var),
                Yn1_view.GetConstColumnView(i_var));
            if (!retval)
            {
              return;
            }
          }
        },
        residual,
        Yn1)(residual, Yn1);
    return retval;
  }
}  // namespace micm
