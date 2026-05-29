// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{
  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
  inline SolverResult AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>::Solve(
      double time_step,
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

    std::size_t max_iter = parameters.max_number_of_steps_;
    const auto TIME_STEP_REDUCTIONS = parameters.time_step_reductions_;

    double h = parameters.h_start_ == 0.0 ? time_step : parameters.h_start_;
    double present_time = 0.0;
    std::size_t n_successful_integrations = 0;
    std::size_t n_convergence_failures = 0;

    auto derived_class_temporary_variables =
        static_cast<BackwardEulerTemporaryVariables<DenseMatrixPolicy>*>(state.temporary_variables_.get());
    auto& yn = derived_class_temporary_variables->Yn_;
    auto& yn1 = state.variables_;  // Yn1 will hold the new solution at the end of the solve
    auto& forcing = derived_class_temporary_variables->forcing_;

    // Ensure Yn starts with the same values as the state variables
    yn.Copy(yn1);

    while (present_time < time_step)
    {
      result.state_ = SolverState::RUNNING;
      bool converged = false;
      std::size_t iterations = 0;

      do
      {
        result.stats_.number_of_steps_++;
        // the first time Yn1 is equal to Yn
        // after the first iteration Yn1 is updated to the new solution
        // so we can use Yn1 to calculate the forcing and jacobian
        // calculate forcing
        forcing.Fill(0.0);
        rates_.AddForcingTerms(state, yn1, forcing);
        result.stats_.function_calls_++;

        // calculate the negative jacobian
        state.jacobian_.Fill(0.0);
        rates_.SubtractJacobianTerms(state, yn1, state.jacobian_);
        result.stats_.jacobian_updates_++;

        // add the inverse of the time step from the diagonal
        state.jacobian_.AddToDiagonal(1 / h);

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
        forcing.ForEach([&](double& f, const double& yn1, const double& yn) { f -= (yn1 - yn) / h; }, yn1, yn);

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
        yn1.ForEach([&](double& yn1, const double& f) { yn1 = std::max(0.0, yn1 + f); }, forcing);

        // if this is the first iteration, we don't need to check for convergence
        if (iterations++ == 0)
        {
          continue;
        }

        // check for convergence
        converged = IsConverged(parameters, forcing, yn1, state.absolute_tolerance_, state.relative_tolerance_);
      } while (!converged && iterations < max_iter);

      if (!converged)
      {
        result.stats_.rejected_++;
        n_successful_integrations = 0;

        if (n_convergence_failures >= TIME_STEP_REDUCTIONS.size())
        {
          present_time += h;
          result.state_ = SolverState::ACCEPTING_UNCONVERGED_INTEGRATION;
          break;
        }
        else
        {
          // if we fail, we need to reset the solution to the last known good solution
          yn1.Copy(yn);
          h *= TIME_STEP_REDUCTIONS[n_convergence_failures++];
        }
      }
      else
      {
        result.stats_.accepted_++;
        result.state_ = SolverState::CONVERGED;
        present_time += h;
        yn.Copy(yn1);

        // when we accept two solutions in a row, we can increase the time step
        n_successful_integrations++;
        if (n_successful_integrations >= 2)
        {
          n_successful_integrations = 0;
          h *= 2.0;
        }
      }
      // Don't let h go past the time step
      h = std::min(h, time_step - present_time);
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
      const std::vector<double>& absolute_tolerance,
      double relative_tolerance)
    requires(!VectorizableDense<DenseMatrixPolicy>)
  {
    double small = parameters.small_;
    double rel_tol = relative_tolerance;
    auto& abs_tol = absolute_tolerance;
    auto residual_iter = residual.AsVector().begin();
    auto yn1_iter = Yn1.AsVector().begin();
    const std::size_t N_ELEM = residual.NumRows() * residual.NumColumns();
    const std::size_t N_VARS = abs_tol.size();
    for (std::size_t i = 0; i < N_ELEM; ++i)
    {
      if (std::abs(*residual_iter) > small && std::abs(*residual_iter) > abs_tol[i % N_VARS] &&
          std::abs(*residual_iter) > rel_tol * std::abs(*yn1_iter))
      {
        return false;
      }
      ++residual_iter, ++yn1_iter;
    }
    return true;
  }

  template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy>
  template<class DenseMatrixPolicy>
  inline bool AbstractBackwardEuler<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>::IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1,
      const std::vector<double>& absolute_tolerance,
      double relative_tolerance)
    requires(VectorizableDense<DenseMatrixPolicy>)
  {
    double small = parameters.small_;
    double rel_tol = relative_tolerance;
    auto& abs_tol = absolute_tolerance;
    auto residual_iter = residual.AsVector().begin();
    auto yn1_iter = Yn1.AsVector().begin();
    const std::size_t N_ELEM = residual.NumRows() * residual.NumColumns();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    const std::size_t N_VARS = abs_tol.size();
    const std::size_t WHOLE_BLOCKS = std::floor(residual.NumRows() / L) * residual.GroupSize();
    // evaluate the rows that fit exactly into the vectorizable dimension (L)
    for (std::size_t i = 0; i < WHOLE_BLOCKS; ++i)
    {
      if (std::abs(*residual_iter) > small && std::abs(*residual_iter) > abs_tol[(i / L) % N_VARS] &&
          std::abs(*residual_iter) > rel_tol * std::abs(*yn1_iter))
      {
        return false;
      }
      ++residual_iter, ++yn1_iter;
    }

    // evaluate the remaining rows
    const std::size_t REMAINING_ROWS = residual.NumRows() % L;
    if (REMAINING_ROWS > 0)
    {
      for (std::size_t y = 0; y < residual.NumColumns(); ++y)
      {
        const std::size_t OFFSET = y * L;
        for (std::size_t i = OFFSET; i < OFFSET + REMAINING_ROWS; ++i)
        {
          if (std::abs(residual_iter[i]) > small && std::abs(residual_iter[i]) > abs_tol[y] &&
              std::abs(residual_iter[i]) > rel_tol * std::abs(yn1_iter[i]))
          {
            return false;
          }
        }
      }
    }
    return true;
  }
}  // namespace micm
