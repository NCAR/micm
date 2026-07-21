// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Smoke test for the RODAS4P tableau
// (RosenbrockSolverParameters::Rodas4PDifferentialAlgebraicRosenbrockParameters):
// on a slaved-equilibrium DAE it must reproduce the analytic solution and agree
// with the standard six-stage RODAS4 method. The tableau's design property —
// avoiding Prothero-Robinson order reduction — is measured by the
// dae_convergence_experiments benchmark harness.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace micm;

namespace
{
  // A -> P (k = 1), slaved equilibrium B = K_eq * A (B algebraic).
  double RunSlavedDae(RosenbrockSolverParameters options, double total_time)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto P = Species("P");
    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, P } };
    const double K_eq = 1.0e-5;

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { P, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0 })
                      .SetPhase(gas_phase)
                      .Build();
    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq",
        B,
        std::vector<StoichSpecies>{ { A, 1.0 } },
        std::vector<StoichSpecies>{ { B, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = 0.0 }));

    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-8);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-14));
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = K_eq;
    state.variables_[0][state.variable_map_.at("P")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    double advanced = 0.0;
    int guard = 0;
    while (advanced < total_time && guard++ < 10000)
    {
      auto result = solver.Solve(total_time - advanced, state);
      EXPECT_EQ(result.state_, SolverState::Converged);
      if (result.state_ != SolverState::Converged)
        break;
      advanced += result.stats_.final_time_;
    }
    // B must stay slaved to A.
    const double a = state.variables_[0][state.variable_map_.at("A")];
    const double b = state.variables_[0][state.variable_map_.at("B")];
    EXPECT_NEAR(b / (K_eq * a), 1.0, 1.0e-6);
    return a;
  }
}  // namespace

TEST(Rodas4P, MatchesAnalyticSolutionAndRodas4)
{
  constexpr double kTotalTime = 5.0;
  const double exact = std::exp(-kTotalTime);

  const double rodas4 = RunSlavedDae(RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters(), kTotalTime);
  const double rodas4p = RunSlavedDae(RosenbrockSolverParameters::Rodas4PDifferentialAlgebraicRosenbrockParameters(), kTotalTime);

  EXPECT_NEAR(rodas4p, exact, 1.0e-6 * exact);
  EXPECT_NEAR(rodas4p, rodas4, 1.0e-6 * exact);
}
