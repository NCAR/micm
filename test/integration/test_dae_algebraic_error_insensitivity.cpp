// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests verifying that the Rosenbrock DAE step-change error estimate for
// algebraic variables makes the solver sensitive to algebraic tolerances
// and prevents overshoot.
//
// Root cause: The embedded error formula Yerror = sum(e_i * K_i) produces
// near-zero entries for algebraic variables because the mass matrix diagonal
// M_ii = 0 zeroes out the inter-stage coupling terms (c/H) * M_ii * K[j].
// Fix: replace Yerror[a] with Ynew[a] - Y[a] for algebraic variables.
//
// System (4 species, 2 equilibria, 1 reaction, 1 conservation):
//   a_gas <-> a_aq   (equilibrium, K1)    -- a_aq is algebraic
//   a_aq  <-> b_aq   (equilibrium, K2)    -- b_aq is algebraic
//   b_aq  -> p       (kinetics, rate k)   -- p is the only differential variable
//   Conservation: a_aq + b_aq + p + a_gas = c_total   -- a_gas is the algebraic balance variable
//
// Analytical solution:
//   a_gas(t) = (C - p(t)) / (1 + K1 + K1*K2)
//   p(t)     = C * (1 - exp(-r*t))   where r = k*K1*K2 / (1 + K1 + K1*K2)
//
// As p approaches c_total, a_gas approaches zero. If the solver overshoots
// p > c_total, the conservation constraint forces a_gas negative.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace micm;

namespace
{
  struct SolveResult
  {
    double A_gas_final_;
    double P_final_;
    double min_A_gas_;
    uint64_t accepted_;
    uint64_t rejected_;
  };

  /// @brief Run a cascade equilibrium system with the given tolerance for the balance variable (a_gas)
  /// @param balance_atol Absolute tolerance for a_gas (the algebraic balance variable)
  /// @param k Rate constant for b_aq -> p reaction
  /// @param K1 Equilibrium constant for a_gas <-> a_aq
  /// @param K2 Equilibrium constant for a_aq <-> b_aq
  SolveResult RunCascadeSystem(double balance_atol, double k, double K1, double K2)
  {
    auto a_gas = Species("A_gas");
    auto a_aq = Species("A_aq");
    auto b_aq = Species("B_aq");
    auto p = Species("P");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a_gas, a_aq, b_aq, p } };

    double c_total = 1.0e-6;

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ b_aq })
                      .SetProducts({ { p, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq1",
        a_aq,
        std::vector<StoichSpecies>{ { a_gas, 1.0 } },
        std::vector<StoichSpecies>{ { a_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K1, .delta_H_ = 0.0 }));
    constraints.push_back(EquilibriumConstraint(
        "eq2",
        b_aq,
        std::vector<StoichSpecies>{ { a_aq, 1.0 } },
        std::vector<StoichSpecies>{ { b_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K2, .delta_H_ = 0.0 }));
    // a_gas is explicitly set as the algebraic balance variable
    constraints.push_back(
        LinearConstraint("mass", a_gas, { { a_aq, 1.0 }, { b_aq, 1.0 }, { p, 1.0 }, { a_gas, 1.0 } }, c_total));

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

    // balance_atol for a_gas (algebraic balance variable under test),
    // moderate atol for equilibrium algebraic variables (a_aq, b_aq),
    // tight atol for the only differential variable (p)
    std::vector<double> atols(4, 1.0e-12);
    atols[gi] = balance_atol;
    atols[ai] = 1.0e-8;
    atols[bi] = 1.0e-8;
    state.SetAbsoluteTolerances(atols);

    // Initial: equilibrium satisfied, p = 0
    double denom = 1.0 + K1 + K1 * K2;
    double ag0 = c_total / denom;
    state.variables_[0][gi] = ag0;
    state.variables_[0][ai] = K1 * ag0;
    state.variables_[0][bi] = K2 * K1 * ag0;
    state.variables_[0][pi] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;

    solver.UpdateStateParameters(state);

    double dt = 30.0;
    double advanced = 0.0;
    double min_a_gas = 1e99;
    uint64_t total_accepted = 0;
    uint64_t total_rejected = 0;

    while (advanced < dt)
    {
      auto result = solver.Solve(dt - advanced, state);
      EXPECT_EQ(result.state_, SolverState::CONVERGED);
      if (result.state_ != SolverState::CONVERGED)
      {
        break;
      }
      advanced += result.stats_.final_time_;
      total_accepted += result.stats_.accepted_;
      total_rejected += result.stats_.rejected_;
      min_a_gas = std::min(min_a_gas, state.variables_[0][gi]);
    }

    return { state.variables_[0][gi], state.variables_[0][pi], min_a_gas, total_accepted, total_rejected };
  }
}  // namespace

/// @brief Prove that the step-change error estimate makes the solver sensitive to algebraic atol.
///
/// With the step-change error injection, tighter atol for the algebraic balance
/// variable (a_gas) should produce more internal steps during the transient.
/// Uses moderate stiffness so the transient is long enough to require many steps.
TEST(DAEAlgebraicError, ErrorSensitiveToBalanceAtol)
{
  // Moderate stiffness: k=100, K1=2, K2=2
  // Time constant ≈ 1/(k*K1*K2/(1+K1+K1*K2)) = 7/(100*4) ≈ 0.018 s
  // A_gas0 = C/7 ≈ 1.43e-7
  auto r_loose = RunCascadeSystem(1e-3, 100.0, 2.0, 2.0);
  auto r_tight = RunCascadeSystem(1e-8, 100.0, 2.0, 2.0);

  // With the fix, tight balance atol should produce more steps
  EXPECT_GT(r_tight.accepted_, r_loose.accepted_) << "Tight atol should require more steps than loose atol. "
                                                << "loose=" << r_loose.accepted_ << " tight=" << r_tight.accepted_;
}

/// @brief Verify that the algebraic balance variable does not go deeply negative.
///
/// Uses an ultra-stiff system (k=1e4, K1=100, K2=50) where the embedded error
/// estimate produces near-zero Yerror for all variables. Without the step-change
/// fix, the solver would accept a huge first step, overshooting p > c_total and
/// forcing a_gas deeply negative.
TEST(DAEAlgebraicError, AlgebraicVariableDoesNotOvershootDeeply)
{
  // Ultra-stiff: transient is ~5e-3 s, solver converges quickly
  auto r = RunCascadeSystem(1e-8, 1.0e4, 100.0, 50.0);

  // a_gas should stay non-negative (or very close to zero).
  // The analytical solution has a_gas >= 0 at all times.
  EXPECT_GE(r.min_A_gas_, -1.0e-8) << "A_gas overshot deeply negative: min_A_gas=" << r.min_A_gas_;

  // After 30s the system should be fully converted: p ≈ c_total, a_gas ≈ 0
  EXPECT_NEAR(r.P_final_, 1.0e-6, 1.0e-8);
  EXPECT_NEAR(r.A_gas_final_, 0.0, 1.0e-10);
}
