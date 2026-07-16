// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for the DAE Rosenbrock algebraic-variable error treatment.
//
// These pin two requirements for the embedded-error treatment now used for
// algebraic rows. See
// docs/superpowers/notes/2026-05-28-dae-algebraic-error-investigation.md.
//
//   FEASIBILITY: an algebraic conservation/balance variable must not be driven
//   into an unphysical (negative) state under stiff conditions. This is the
//   property that motivated the former step-change override.
//
//   EFFICIENCY: a *slaved* algebraic variable that legitimately spans orders of
//   magnitude must not inflate the solver step count when its absolute tolerance
//   is tightened. The former O(H) raw-step-change override made this cost grow
//   dramatically; the embedded LTE must keep the step count insensitive to a
//   tolerance on a slaved algebraic variable.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace micm;

namespace
{
  struct CascadeResult
  {
    double final_P;
    double min_balance;  // smallest value the algebraic balance variable A_gas reached
    std::uint64_t accepted;
    bool converged;
  };

  // Cascade equilibrium + conservation system (same family as
  // test_dae_algebraic_error_insensitivity.cpp):
  //   A_gas <-> A_aq  (equilibrium K1, A_aq algebraic)
  //   A_aq  <-> B_aq  (equilibrium K2, B_aq algebraic)
  //   B_aq  -> P      (kinetics rate k, P differential)
  //   conservation: A_gas + A_aq + B_aq + P = C_total  (A_gas algebraic balance)
  // As P -> C_total the balance variable A_gas -> 0; overshoot would drive it negative.
  CascadeResult RunCascade(double balance_atol, double k, double K1, double K2)
  {
    auto A_gas = Species("A_gas");
    auto A_aq = Species("A_aq");
    auto B_aq = Species("B_aq");
    auto P = Species("P");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A_gas, A_aq, B_aq, P } };

    double C_total = 1.0e-6;

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ B_aq })
                      .SetProducts({ { P, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq1",
        A_aq,
        std::vector<StoichSpecies>{ { A_gas, 1.0 } },
        std::vector<StoichSpecies>{ { A_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K1, .delta_H_ = 0.0 }));
    constraints.push_back(EquilibriumConstraint(
        "eq2",
        B_aq,
        std::vector<StoichSpecies>{ { A_aq, 1.0 } },
        std::vector<StoichSpecies>{ { B_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K2, .delta_H_ = 0.0 }));
    constraints.push_back(
        LinearConstraint("mass", A_gas, { { A_aq, 1.0 }, { B_aq, 1.0 }, { P, 1.0 }, { A_gas, 1.0 } }, C_total));

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);

    auto& vm = state.variable_map_;
    std::size_t gi = vm.at("A_gas");
    std::size_t ai = vm.at("A_aq");
    std::size_t bi = vm.at("B_aq");
    std::size_t pi = vm.at("P");

    // Tight atol for the differential variable (P); the balance variable (A_gas)
    // tolerance is the knob under test; moderate for the other algebraic species.
    std::vector<double> atols(4, 1.0e-12);
    atols[gi] = balance_atol;
    atols[ai] = 1.0e-8;
    atols[bi] = 1.0e-8;
    state.SetAbsoluteTolerances(atols);

    double denom = 1.0 + K1 + K1 * K2;
    double Ag0 = C_total / denom;
    state.variables_[0][gi] = Ag0;
    state.variables_[0][ai] = K1 * Ag0;
    state.variables_[0][bi] = K2 * K1 * Ag0;
    state.variables_[0][pi] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    CascadeResult out{ 0.0, 1e99, 0, true };
    double dt = 30.0;
    double advanced = 0.0;
    int guard = 0;
    while (advanced < dt && guard++ < 100000)
    {
      auto result = solver.Solve(dt - advanced, state);
      if (result.state_ != SolverState::Converged)
      {
        out.converged = false;
        break;
      }
      advanced += result.stats_.final_time_;
      out.accepted += result.stats_.accepted_;
      out.min_balance = std::min(out.min_balance, state.variables_[0][gi]);
    }
    out.final_P = state.variables_[0][pi];
    return out;
  }

  struct DecayResult
  {
    std::uint64_t accepted;
    bool converged;
  };

  // Slow decay coupled to a slaved equilibrium variable, integrated over many
  // decades (Robertson-like): this is the regime where the former step-change
  // error estimate inflated step counts, because the algebraic variable
  // B = K_eq*A slowly sweeps orders of magnitude as A decays.
  //   A -> P            (rate k, A and P differential; A+P conserved)
  //   K_eq*A - B = 0    (B algebraic, slaved to A; B0 = K_eq ~ 1e-5 scale)
  // Tightening atol_B throttles H to keep dB ~ atol_B per step over the whole
  // decay, inflating the accepted-step count for no accuracy gain.
  DecayResult RunDecayEquilibrium(double balance_atol)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto P = Species("P");
    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, P } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { P, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    double K_eq = 1.0e-5;
    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "slaved",
        B,
        std::vector<StoichSpecies>{ { A, 1.0 } },
        std::vector<StoichSpecies>{ { B, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = 0.0 }));

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    auto& vm = state.variable_map_;
    std::vector<double> atols(3, 1.0e-12);
    atols[vm.at("B")] = balance_atol;
    state.SetAbsoluteTolerances(atols);

    state.variables_[0][vm.at("A")] = 1.0;
    state.variables_[0][vm.at("B")] = K_eq * 1.0;
    state.variables_[0][vm.at("P")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    DecayResult out{ 0, true };
    double dt = 50.0;
    double advanced = 0.0;
    int guard = 0;
    while (advanced < dt && guard++ < 100000)
    {
      auto result = solver.Solve(dt - advanced, state);
      if (result.state_ != SolverState::Converged && result.state_ != SolverState::ConvergenceExceededMaxSteps)
      {
        out.converged = false;
        break;
      }
      out.accepted += result.stats_.accepted_;
      if (result.stats_.final_time_ <= 0.0)
      {
        out.converged = false;
        break;
      }
      advanced += result.stats_.final_time_;
    }
    return out;
  }
}  // namespace

// FEASIBILITY CONTRACT (active): under aggressive stiffness, the algebraic
// balance variable must not be driven negative. Any future change to the
// algebraic error handling must keep this passing.
TEST(DAEAlgebraicStepEconomy, BalanceVariableStaysNonNegativeUnderStiffness)
{
  // Ultra-stiff cascade: fast drain to P, strong equilibria.
  auto r = RunCascade(/*balance_atol=*/1.0e-8, /*k=*/1.0e4, /*K1=*/100.0, /*K2=*/50.0);
  ASSERT_TRUE(r.converged) << "solver failed to converge";
  EXPECT_GE(r.min_balance, -1.0e-8) << "algebraic balance variable A_gas overshot negative: min=" << r.min_balance;
  // Fully converted at 30 s: P ~ C_total, balance ~ 0.
  EXPECT_NEAR(r.final_P, 1.0e-6, 1.0e-8);
}

// EFFICIENCY CONTRACT: tightening the balance variable's absolute tolerance
// must NOT substantially inflate the accepted step count. The balance variable
// is slaved by the conservation constraint; its accuracy is determined by the
// differential variable and the constraint, not by its own atol. The solver
// uses the method's embedded local truncation error for algebraic rows (which
// is ~0 for slaved variables), so a tighter algebraic atol does not throttle the
// step size. (Before the algebraic error handling was corrected, the O(H)
// step-change estimate made this fail ~100x.)
TEST(DAEAlgebraicStepEconomy, TightBalanceAtolDoesNotInflateSteps)
{
  auto loose = RunDecayEquilibrium(/*balance_atol=*/1.0e-2);
  auto tight = RunDecayEquilibrium(/*balance_atol=*/1.0e-10);
  ASSERT_TRUE(loose.converged && tight.converged);

  // The algebraic variable B = K_eq*A is slaved; its accuracy is determined by A
  // and the constraint, not by its own atol. A correct algebraic-error treatment
  // should not pay a large step-count penalty for a tighter atol on it. Under the
  // former step-change estimate, tight atol inflated the count many-fold.
  EXPECT_LE(static_cast<double>(tight.accepted), 2.0 * static_cast<double>(loose.accepted))
      << "tightening balance atol inflated steps: loose=" << loose.accepted << " tight=" << tight.accepted;
}
