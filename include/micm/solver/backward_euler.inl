// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
enum class MicmBackwardEulerErrc
{
  FailedToConverge = 1
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmBackwardEulerErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class BackwardEulerErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return "MICM BackwardEuler";
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmBackwardEulerErrc>(ev))
      {
        case MicmBackwardEulerErrc::FailedToConverge: return "Failed to converge";
        default: return "Unknown error";
      }
    }
  };

  const BackwardEulerErrorCategory backwardEulerErrorCategory{};
}  // namespace

inline std::error_code make_error_code(MicmBackwardEulerErrc e)
{
  return { static_cast<int>(e), backwardEulerErrorCategory };
}

namespace micm
{
  template<class RatesPolicy, class LinearSolverPolicy>
  inline SolverResult BackwardEuler<RatesPolicy, LinearSolverPolicy>::Solve(double time_step, auto& state) const
  {
    // A fully implicit euler implementation is given by the following equation:
    // y_{n+1} = y_n + H * f(t_{n+1}, y_{n+1})
    // This is a root finding problem because you need to know y_{n+1} to compute f(t_{n+1}, y_{n+1})
    // you need to solve the equation y_{n+1} - y_n - H f(t_{n+1}, y_{n+1}) = 0
    // A series of time step reductions are used after failed solves to try to find a solution
    // These reductions are controlled by the time_step_reductions parameter in the solver parameters
    // if the last attempt to reduce the timestep fails,
    // accept the current H but do not update the Yn vector

    using MatrixPolicy = decltype(state.variables_);

    SolverResult result;

    std::size_t max_iter = parameters_.max_number_of_steps_;
    const auto time_step_reductions = parameters_.time_step_reductions_;

    double H = time_step;
    double t = 0.0;
    std::size_t n_successful_integrations = 0;
    std::size_t n_convergence_failures = 0;

    auto derived_class_temporary_variables =
        static_cast<BackwardEulerTemporaryVariables<MatrixPolicy>*>(state.temporary_variables_.get());
    auto& Yn = derived_class_temporary_variables->Yn_;
    auto& Yn1 = state.variables_;  // Yn1 will hold the new solution at the end of the solve
    auto& forcing = derived_class_temporary_variables->forcing_;

    while (t < time_step)
    {
      result.state_ = SolverState::Running;
      bool converged = false;
      std::size_t iterations = 0;

      if (result.stats_.number_of_steps_ == 0)
      {
        Yn.Copy(Yn1);
      }
      else
      {
        Yn1.Copy(Yn);
      }

      do
      {
        result.stats_.number_of_steps_++;
        // the first time Yn1 is equal to Yn
        // after the first iteration Yn1 is updated to the new solution
        // so we can use Yn1 to calculate the forcing and jacobian
        // calculate forcing
        forcing.Fill(0.0);
        rates_.AddForcingTerms(state.rate_constants_, Yn1, forcing);
        result.stats_.function_calls_++;

        // calculate the negative jacobian
        state.jacobian_.Fill(0.0);
        rates_.SubtractJacobianTerms(state.rate_constants_, Yn1, state.jacobian_);
        result.stats_.jacobian_updates_++;

        // add the inverse of the time step from the diagonal
        state.jacobian_.AddToDiagonal(1 / H);

        // We want to solve this equation for a zero
        // (y_{n+1} - y_n) / H = f(t_{n+1}, y_{n+1})

        // try to find the root by factoring and solving the linear system
        linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
        result.stats_.decompositions_++;

        // forcing_blk in camchem
        // residual = forcing - (Yn1 - Yn) / H
        // since forcing is only used once, we can reuse it to store the residual
        forcing.ForEach([&](double& f, const double& yn1, const double& yn) { f -= (yn1 - yn) / H; }, Yn1, Yn);

        // the result of the linear solver will be stored in forcing
        // this represents the change in the solution
        linear_solver_.Solve(forcing, state.lower_matrix_, state.upper_matrix_);
        result.stats_.solves_++;

        // solution_blk in camchem
        // Yn1 = Yn1 + residual;
        // always make sure the solution is positive regardless of which iteration we are on
        Yn1.ForEach([&](double& yn1, const double& f) { yn1 = std::max(0.0, yn1 + f); }, forcing);

        // if this is the first iteration, we don't need to check for convergence
        if (iterations++ == 0)
          continue;

        // check for convergence
        converged = IsConverged(parameters_, forcing, Yn1, state.absolute_tolerance_, state.relative_tolerance_);
      } while (!converged && iterations < max_iter);

      if (!converged)
      {
        result.stats_.rejected_++;
        n_successful_integrations = 0;

        if (n_convergence_failures >= time_step_reductions.size())
        {
          t += H;
          result.state_ = SolverState::AcceptingUnconvergedIntegration;
          break;
        }
        else
        {
          H *= time_step_reductions[n_convergence_failures++];
        }
      }
      else
      {
        result.stats_.accepted_++;
        result.state_ = SolverState::Converged;
        t += H;
        Yn = Yn1;

        // when we accept two solutions in a row, we can increase the time step
        n_successful_integrations++;
        if (n_successful_integrations >= 2)
        {
          n_successful_integrations = 0;
          H *= 2.0;
        }
      }
      // Don't let H go past the time step
      H = std::min(H, time_step - t);
    }

    result.final_time_ = t;
    return result;
  }

  template<class RatesPolicy, class LinearSolverPolicy>
  template<class DenseMatrixPolicy>
  inline bool BackwardEuler<RatesPolicy, LinearSolverPolicy>::IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1, std::vector<double>& absolute_tolerance, double relative_tolerance) requires(!VectorizableDense<DenseMatrixPolicy>)
  {
    double small = parameters.small_;
    double rel_tol = relative_tolerance;
    auto& abs_tol = absolute_tolerance;
    auto residual_iter = residual.AsVector().begin();
    auto Yn1_iter = Yn1.AsVector().begin();
    const std::size_t n_elem = residual.NumRows() * residual.NumColumns();
    const std::size_t n_vars = abs_tol.size();
    for (std::size_t i = 0; i < n_elem; ++i)
    {
      if (std::abs(*residual_iter) > small && std::abs(*residual_iter) > abs_tol[i % n_vars] &&
          std::abs(*residual_iter) > rel_tol * std::abs(*Yn1_iter))
      {
        return false;
      }
      ++residual_iter, ++Yn1_iter;
    }
    return true;
  }

  template<class RatesPolicy, class LinearSolverPolicy>
  template<class DenseMatrixPolicy>
  inline bool BackwardEuler<RatesPolicy, LinearSolverPolicy>::IsConverged(
      const BackwardEulerSolverParameters& parameters,
      const DenseMatrixPolicy& residual,
      const DenseMatrixPolicy& Yn1, std::vector<double>& absolute_tolerance, double relative_tolerance) requires(VectorizableDense<DenseMatrixPolicy>)
  {
    double small = parameters.small_;
    double rel_tol = relative_tolerance;
    auto& abs_tol = absolute_tolerance;
    auto residual_iter = residual.AsVector().begin();
    auto Yn1_iter = Yn1.AsVector().begin();
    const std::size_t n_elem = residual.NumRows() * residual.NumColumns();
    const std::size_t L = residual.GroupVectorSize();
    const std::size_t n_vars = abs_tol.size();
    const std::size_t whole_blocks = std::floor(residual.NumRows() / L) * residual.GroupSize();
    // evaluate the rows that fit exactly into the vectorizable dimension (L)
    for (std::size_t i = 0; i < whole_blocks; ++i)
    {
      if (std::abs(*residual_iter) > small && std::abs(*residual_iter) > abs_tol[(i / L) % n_vars] &&
          std::abs(*residual_iter) > rel_tol * std::abs(*Yn1_iter))
      {
        return false;
      }
      ++residual_iter, ++Yn1_iter;
    }

    // evaluate the remaining rows
    const std::size_t remaining_rows = residual.NumRows() % L;
    if (remaining_rows > 0)
    {
      for (std::size_t y = 0; y < residual.NumColumns(); ++y)
      {
        const std::size_t offset = y * L;
        for (std::size_t i = offset; i < offset + remaining_rows; ++i)
        {
          if (std::abs(residual_iter[i]) > small && std::abs(residual_iter[i]) > abs_tol[y] &&
              std::abs(residual_iter[i]) > rel_tol * std::abs(Yn1_iter[i]))
          {
            return false;
          }
        }
      }
    }
    return true;
  }
}  // namespace micm
