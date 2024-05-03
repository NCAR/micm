/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
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
        case MicmBackwardEulerErrc::FailedToConverge:
          return "Failed to converge";
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
  inline BackwardEuler::BackwardEuler()
  {
  }

  inline void BackwardEuler::Solve(double time_step, auto& state, auto linear_solver, auto process_set, const std::vector<micm::Process>& processes, auto jacobian_diagonal_elements)
  {
    // A fully implicit euler implementation is given by the following equation:
    // y_{n+1} = y_n + H * f(t_{n+1}, y_{n+1})
    // This is a root finding problem because you need to know y_{n+1} to compute f(t_{n+1}, y_{n+1})
    // you need to solve the equation y_{n+1} - y_n - H f(t_{n+1}, y_{n+1}) = 0
    // We will also use the same logic used by cam-chem to determine the time step
    // That scheme is this: 
    // Start with H = time_step
    // if that fails, try H = H/2
    // if that fails, try H = H/10
    // if that fails, return whatever the current integration with the current H is

    double H = time_step;
    double t = 0.0;
    double tol = 1e-4;
    double small = 1.0e-40;
    auto temp = state.variables_;
    auto forcing = state.variables_;
    auto residual = state.variables_;
    bool singular = false;
    std::size_t max_iter = 11;
    std::size_t n_successful_integrations = 0;
    std::size_t n_convergence_failures = 0;
    std::size_t cut_in_half_limit = 4;
    std::size_t cut_in_tenth_limit = 1;
    const std::array<double, 5> time_step_reductions{ 0.5, 0.5, 0.5, 0.5, 0.1 };

    auto Yn = state.variables_;
    auto Yn1 = state.variables_;

    Process::UpdateState(processes, state);

    while(t < time_step) {
      bool converged = true;
      std::size_t iterations = 0;

      do {
        // calculate forcing
        std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
        process_set.AddForcingTerms(state.rate_constants_, Yn, forcing);

        // calculate jacobian
        std::fill(state.jacobian_.AsVector().begin(), state.jacobian_.AsVector().end(), 0.0);
        process_set.SubtractJacobianTerms(state.rate_constants_, Yn, state.jacobian_);

        // subtract the inverse of the time step from the diagonal
        // TODO: handle vectorized jacobian matrix
        for(auto& jac : state.jacobian_.AsVector()) {
          jac *= -1;
        }
        for (std::size_t i_block = 0; i_block < state.jacobian_.Size(); ++i_block)
        {
          auto jacobian_vector = std::next(state.jacobian_.AsVector().begin(), i_block * state.jacobian_.FlatBlockSize());
          for (const auto& i_elem : jacobian_diagonal_elements)
            jacobian_vector[i_elem] -= 1 / H;
        }

        // We want to solve this equation for a zero
        // (y_{n+1} - y_n) / H = f(t_{n+1}, y_{n+1})

        // try to find the root by factoring and solving the linear system
        linear_solver.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_, singular);

        auto yn1_iter = Yn1.begin();
        auto yn_iter = Yn.begin();
        auto residual_iter = residual.begin();
        auto forcing_iter = forcing.begin();
        // forcing_blk in camchem
        // residual = (Yn1 - Yn) / H - forcing;
        for(; yn1_iter != Yn1.end(); ++yn1_iter, ++yn_iter, ++residual_iter, ++forcing_iter) {
          *residual_iter = (*yn1_iter - *yn_iter) / H - *forcing_iter;
        }

        // the result of the linear solver will be stored in temp
        // this represnts the change in the solution
        linear_solver.Solve(residual, temp, state.lower_matrix_, state.upper_matrix_);

        // solution_blk in camchem
        // Yn1 = Yn1 + temp;
        auto temp_iter = temp.begin();
        yn1_iter = Yn1.begin();

        for(; temp_iter != temp.end(); ++temp_iter, ++yn1_iter) {
          *yn1_iter += *temp_iter;
        }

        // always make sure the solution is positive regardless of which iteration we are on
        // remove negatives, accept this solution and continue
        yn1_iter = Yn1.begin();
        yn_iter = Yn.begin();
        for(; yn1_iter != Yn1.end(); ++yn1_iter, ++yn_iter) {
          *yn_iter = *yn1_iter = std::max(0.0, *yn1_iter);
        }

        if (iterations == 0) {
          // we always want to accept the first iteration
          continue;
        }

        // check for convergence
        temp_iter = temp.begin();
        yn1_iter = Yn1.begin();

        converged = true;
        // convergence happens when the absolute value of the change to the solution
        // is less than a tolerance times the absolute value of the solution
        for(; temp_iter != temp.end(); ++temp_iter, ++yn1_iter) {
          // changes that are much smaller than the tolerance are negligible and we assume can be accepted
          if (std::abs(*temp_iter) > small) {
            if (std::abs(*temp_iter) > tol * std::abs(*yn1_iter)) {
              converged = false;
              break;
            }
          }
        }

        ++iterations;
        // correct this check
      } while(!converged && iterations < max_iter);

      if (!converged) {
        n_successful_integrations = 0;

        if (n_convergence_failures >= time_step_reductions.size()) {
          // we have failed to converge, accept the solution
          // TODO: continue on with the current solution to get the full solution
          n_convergence_failures = 0;
          // give_up = true;
          t += H;
          break;
          throw std::system_error(make_error_code(MicmBackwardEulerErrc::FailedToConverge), "Failed to converge");
        };

        H *= time_step_reductions[n_convergence_failures++];
      }
      else {
        t += H;

        // when we accept two solutions in a row, we can increase the time step
        n_successful_integrations++;
        n_convergence_failures = 0;
        if (n_successful_integrations >= 2) {
          H *= 2.0;
        }
      }

      // Don't let H go past the time step
      H = std::min(H, time_step - t);
    }

    state.variables_ = Yn1;
  }
}
