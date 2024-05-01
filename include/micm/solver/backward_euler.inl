/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
enum class MicmBackwardEulerErrc
{
  SubStepTooSmall = 1,
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
        case MicmBackwardEulerErrc::SubStepTooSmall:
          return "Substep too small";
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
      : processes_()
  {
  }

  inline BackwardEuler::BackwardEuler(
        const System& system,
        const std::vector<Process>& processes)
      : processes_()
  {
  }

  enum class LastTimeStepCut
  {
    None,
    Half,
    Tenth
  };

  inline void BackwardEuler::Solve(double time_step, auto state, auto linear_solver, auto process_set)
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
    double tol = 1.0e-6;
    auto y = state.variables_;
    auto temp = y;
    auto forcing = y;
    auto residual = y;
    enum LastTimeStepCut last_time_step_cut = LastTimeStepCut::None;

    auto Yn = y;
    auto Yn1 = y;

    Process::UpdateState(processes_, state);

    bool singular = false;

    while(t < time_step) {

      bool converged = true;

      do {
        // calculate jacobian
        std::fill(state.jacobian_.AsVector().begin(), state.jacobian_.AsVector().end(), 0.0);
        process_set.SubtractJacobianTerms(state.rate_constants_, Yn, state.jacobian_);

        // calculate forcing
        std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
        process_set.AddForcingTerms(state.rate_constants_, Yn, forcing);

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
        // Yn1 = Yn + temp;
        auto temp_iter = temp.begin();
        yn1_iter = Yn1.begin();
        yn_iter = Yn.begin();

        for(; temp_iter != temp.end(); ++temp_iter, ++yn1_iter, ++yn_iter) {
          *yn1_iter = *yn_iter + *temp_iter;
        }

        // check for convergence
        temp_iter = temp.begin();
        yn1_iter = Yn1.begin();

        // convergence happens when the absolute value of the change to the solution
        // is less than a tolerance times the absolute value of the solution
        for(; temp_iter != temp.end(); ++temp_iter, ++yn1_iter) {
          if (std::abs(*temp_iter) > tol * std::abs(*yn1_iter)) {
            converged = false;
          }
        }

        if (!converged) {
          switch (last_time_step_cut) {
            case LastTimeStepCut::None:
              last_time_step_cut = LastTimeStepCut::Half;
              H /= 2.0;
              break;
            case LastTimeStepCut::Half:
              last_time_step_cut = LastTimeStepCut::Tenth;
              H /= 10;
              break;
            case LastTimeStepCut::Tenth:
              throw std::system_error(make_error_code(MicmBackwardEulerErrc::SubStepTooSmall));
          }
        }
        else {
          t += H;
        }
      } while(!converged);
    }

    state.variables_ = Yn1;
  }
}