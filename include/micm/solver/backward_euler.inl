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
  inline SolverResult BackwardEuler<RatesPolicy, LinearSolverPolicy>::Solve(double time_step, auto& state)
  {
    // A fully implicit euler implementation is given by the following equation:
    // y_{n+1} = y_n + H * f(t_{n+1}, y_{n+1})
    // This is a root finding problem because you need to know y_{n+1} to compute f(t_{n+1}, y_{n+1})
    // you need to solve the equation y_{n+1} - y_n - H f(t_{n+1}, y_{n+1}) = 0
    // We will also use the same logic used by cam-chem to determine the time step
    // That scheme is this:
    // Start with H = time_step
    // if that fails, try H = H/2 several times
    // if that fails, try H = H/10 once
    // if that fails, accept the current H but do not update the Yn vector

    // TODO populate the result before returning it
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
        // the first time Yn1 is equal to Yn
        // after the first iteration Yn1 is updated to the new solution
        // so we can use Yn1 to calculate the forcing and jacobian
        // calculate forcing
        std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
        rates_.AddForcingTerms(state.rate_constants_, Yn1, forcing);

        // calculate jacobian
        std::fill(state.jacobian_.AsVector().begin(), state.jacobian_.AsVector().end(), 0.0);
        rates_.SubtractJacobianTerms(state.rate_constants_, Yn1, state.jacobian_);

        // subtract the inverse of the time step from the diagonal
        // TODO: handle vectorized jacobian matrix
        for (auto& jac : state.jacobian_.AsVector())
        {
          jac *= -1;
        }
        for (std::size_t i_block = 0; i_block < state.jacobian_.NumberOfBlocks(); ++i_block)
        {
          auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_block * state.jacobian_.FlatBlockSize());
          for (const auto& i_elem : jacobian_diagonal_elements_)
            jacobian_vector[i_elem] -= 1 / H;
        }

        // We want to solve this equation for a zero
        // (y_{n+1} - y_n) / H = f(t_{n+1}, y_{n+1})

        // try to find the root by factoring and solving the linear system
        linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_, singular);

        auto yn1_iter = Yn1.begin();
        auto yn_iter = Yn.begin();
        auto forcing_iter = forcing.begin();
        // forcing_blk in camchem
        // residual = (Yn1 - Yn) / H - forcing;
        // since forcing is only used once, we can reuse it to store the residual
        for (; yn1_iter != Yn1.end(); ++yn1_iter, ++yn_iter, ++forcing_iter)
        {
          *forcing_iter = (*yn1_iter - *yn_iter) / H - *forcing_iter;
        }

        // the result of the linear solver will be stored in forcing
        // this represents the change in the solution
        linear_solver_.Solve(forcing, forcing, state.lower_matrix_, state.upper_matrix_);

        // solution_blk in camchem
        // Yn1 = Yn1 + residual;
        // always make sure the solution is positive regardless of which iteration we are on
        forcing_iter = forcing.begin();
        yn1_iter = Yn1.begin();
        for (; forcing_iter != forcing.end(); ++forcing_iter, ++yn1_iter)
        {
          *yn1_iter = std::max(0.0, *yn1_iter + *forcing_iter);
        }

        // if this is the first iteration, we don't need to check for convergence
        if (iterations++ == 0)
          continue;

        // check for convergence
        forcing_iter = forcing.begin();
        yn1_iter = Yn1.begin();

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
        n_successful_integrations = 0;

        if (n_convergence_failures >= time_step_reductions.size())
        {
          t += H;
          std::cerr << "Failed to converge too many times in a row. Accepting the current integration and continuing on.\n";
          break;
        }
        else {
          std::cerr << "Failed to converge. Reducing the time step.\n";
          H *= time_step_reductions[n_convergence_failures++];
        }

      }
      else
      {
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
