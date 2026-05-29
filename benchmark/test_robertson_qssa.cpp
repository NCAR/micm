// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "robertson_qssa_constraint.hpp"
#include "robertson_system.hpp"

#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace
{
  // Integrate Robertson from t=0 to t_end, returning {A,B,C} at t_end.
  // If qssa==true, B is algebraic via the QSSA constraint; otherwise full ODE.
  std::vector<double> Integrate(bool qssa, double k1, double k2, double k3, double rtol, double t_end)
  {
    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

    // The constraint is copied by value into AddExternalModel; keep it in scope
    // through the build expression. Build in one unbroken chain (never store the
    // intermediate builder in `auto` — the chained setters return a base
    // SolverBuilder& and would slice). Both branches produce the SAME solver
    // type because AddExternalModel returns SolverBuilder& (type-erased), so the
    // ternary is well-typed.
    robertson::QssaConstraint constraint(k1, k2, k3);

    auto solver =
        qssa ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                   .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                   .SetReactions(sys.processes)
                   .AddExternalModel(constraint)
                   .SetReorderState(false)
                   .Build()
             : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                   .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                   .SetReactions(sys.processes)
                   .SetReorderState(false)
                   .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, rtol * 1e-2));
    state.SetCustomRateParameter("r1", k1);
    state.SetCustomRateParameter("r2", k2);
    state.SetCustomRateParameter("r3", k3);

    auto map = state.variable_map_;
    state.variables_[0][map.at("A")] = 1.0;
    state.variables_[0][map.at("B")] = qssa ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
    state.variables_[0][map.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 272.5;
    state.conditions_[0].pressure_ = 101253.3;
    state.conditions_[0].air_density_ = 1e6;
    solver.UpdateStateParameters(state);

    // Integrate the full interval in a resume loop. A single Solve call over the
    // whole 1e3 s interval can hit the per-call step cap (max_number_of_steps_),
    // returning ConvergenceExceededMaxSteps with final_time_ reflecting the
    // progress made; the loop then re-calls on the remaining interval. This
    // matches the resume pattern in test_analytical_robertson (which does not
    // assert per-call convergence). Treat the step-cap state as a valid resume
    // signal and only fail on hard errors (singular/NaN/Inf/too-small step).
    double done = 0.0;
    while (done < t_end)
    {
      auto result = solver.Solve(t_end - done, state);
      EXPECT_TRUE(
          result.state_ == micm::SolverState::Converged ||
          result.state_ == micm::SolverState::ConvergenceExceededMaxSteps)
          << "unexpected solver state " << static_cast<int>(result.state_);
      EXPECT_GT(result.stats_.final_time_, 0.0) << "no progress made; would loop forever";
      done += result.stats_.final_time_;
    }
    return { state.variables_[0][map.at("A")],
             state.variables_[0][map.at("B")],
             state.variables_[0][map.at("C")] };
  }
}  // namespace

// The QSSA-reduced DAE should match the full ODE for A and C after the transient.
TEST(RobertsonQssa, MatchesFullOdePostTransient)
{
  const double k1 = 0.04, k2 = 3e7, k3 = 1e4, t_end = 1.0e3;
  auto ode = Integrate(false, k1, k2, k3, 1e-10, t_end);  // tight-tol reference
  auto dae = Integrate(true, k1, k2, k3, 1e-6, t_end);

  auto rel = [](double got, double ref) { return std::abs(got - ref) / (std::abs(ref) + 1e-30); };
  EXPECT_LT(rel(dae[0], ode[0]), 1e-3) << "A: dae=" << dae[0] << " ode=" << ode[0];
  EXPECT_LT(rel(dae[2], ode[2]), 1e-3) << "C: dae=" << dae[2] << " ode=" << ode[2];
}

// The QSSA constraint residual G = k1*A - k2*B^2 - k3*B*C must be ~0 at the end.
TEST(RobertsonQssa, ConstraintResidualNearZero)
{
  const double k1 = 0.04, k2 = 3e7, k3 = 1e4;
  auto dae = Integrate(true, k1, k2, k3, 1e-6, 1.0e3);
  double g = k1 * dae[0] - k2 * dae[1] * dae[1] - k3 * dae[1] * dae[2];
  EXPECT_NEAR(g, 0.0, 1e-6) << "residual G=" << g;
}
