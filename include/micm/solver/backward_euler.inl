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

    SolverResult result;

    double small = parameters_.small;
    std::size_t max_iter = parameters_.max_number_of_steps_;
    const auto time_step_reductions = parameters_.time_step_reductions;

    double H = time_step;
    double t = 0.0;
    std::size_t n_successful_integrations = 0;
    std::size_t n_convergence_failures = 0;

    bool singular = false;

    auto Yn = state.variables_;
    auto Yn1 = state.variables_;
    auto forcing = state.variables_;

    while (t < time_step)
    {
      result.state_ = SolverState::Running;
      bool converged = false;
      std::size_t iterations = 0;

      Yn1 = Yn;

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
        linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_, singular);
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
        Yn1.ForEach([&](double& yn1, const double& f) { yn1 = std::max( 0.0, yn1 + f ); }, forcing);

        // if this is the first iteration, we don't need to check for convergence
        if (iterations++ == 0)
          continue;

        // check for convergence
        auto forcing_iter = forcing.begin();
        auto yn1_iter = Yn1.begin();

        // convergence happens when the absolute value of the change to the solution
        // is less than a tolerance times the absolute value of the solution
        auto abs_tol_iter = parameters_.absolute_tolerance_.begin();
        do
        {
          // changes that are much smaller than the tolerance are negligible and we assume can be accepted
          converged = (std::abs(*forcing_iter) <= small) || (std::abs(*forcing_iter) <= *abs_tol_iter) ||
                      (std::abs(*forcing_iter) <= parameters_.relative_tolerance_ * std::abs(*yn1_iter));
          ++forcing_iter, ++yn1_iter, ++abs_tol_iter;
        } while (converged && forcing_iter != forcing.end());
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

    state.variables_ = Yn1;
    result.final_time_ = t;
    return result;
  }

}  // namespace micm
