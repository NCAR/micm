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
//   A_gas <-> A_aq   (equilibrium, K1)    -- A_aq is algebraic
//   A_aq  <-> B_aq   (equilibrium, K2)    -- B_aq is algebraic
//   B_aq  -> P       (kinetics, rate k)   -- P is the only differential variable
//   Conservation: A_aq + B_aq + P + A_gas = C_total   -- A_gas is algebraic (last term)
//
// Analytical solution:
//   A_gas(t) = (C - P(t)) / (1 + K1 + K1*K2)
//   P(t)     = C * (1 - exp(-r*t))   where r = k*K1*K2 / (1 + K1 + K1*K2)
//
// As P approaches C_total, A_gas approaches zero. If the solver overshoots
// P > C_total, the conservation constraint forces A_gas negative.

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
    double A_gas_final;
    double P_final;
    double min_A_gas;
    uint64_t accepted;
    uint64_t rejected;
  };

  /// @brief Run a cascade equilibrium system with the given tolerance for the balance variable (A_gas)
  /// @param balance_atol Absolute tolerance for A_gas (the algebraic balance variable)
  /// @param k Rate constant for B_aq -> P reaction
  /// @param K1 Equilibrium constant for A_gas <-> A_aq
  /// @param K2 Equilibrium constant for A_aq <-> B_aq
  SolveResult RunCascadeSystem(double balance_atol, double k, double K1, double K2)
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
                      .SetRateConstant(ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                      .SetPhase(gas_phase)
                      .Build();

    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq1",
        std::vector<StoichSpecies>{ { A_gas, 1.0 } },
        std::vector<StoichSpecies>{ { A_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref = K1, .delta_H = 0.0 }));
    constraints.push_back(EquilibriumConstraint(
        "eq2",
        std::vector<StoichSpecies>{ { A_aq, 1.0 } },
        std::vector<StoichSpecies>{ { B_aq, 1.0 } },
        VantHoffParam{ .K_HLC_ref = K2, .delta_H = 0.0 }));
    // A_gas is LAST so it becomes the algebraic balance variable
    constraints.push_back(LinearConstraint("mass", { { A_aq, 1.0 }, { B_aq, 1.0 }, { P, 1.0 }, { A_gas, 1.0 } }, C_total));

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
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

    // balance_atol for A_gas (algebraic balance variable under test),
    // moderate atol for equilibrium algebraic variables (A_aq, B_aq),
    // tight atol for the only differential variable (P)
    std::vector<double> atols(4, 1.0e-12);
    atols[gi] = balance_atol;
    atols[ai] = 1.0e-8;
    atols[bi] = 1.0e-8;
    state.SetAbsoluteTolerances(atols);

    // Initial: equilibrium satisfied, P = 0
    double denom = 1.0 + K1 + K1 * K2;
    double Ag0 = C_total / denom;
    state.variables_[0][gi] = Ag0;
    state.variables_[0][ai] = K1 * Ag0;
    state.variables_[0][bi] = K2 * K1 * Ag0;
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
      EXPECT_EQ(result.state_, SolverState::Converged);
      if (result.state_ != SolverState::Converged)
        break;
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
/// variable (A_gas) should produce more internal steps during the transient.
/// Uses moderate stiffness so the transient is long enough to require many steps.
TEST(DAEAlgebraicError, ErrorSensitiveToBalanceAtol)
{
  // Moderate stiffness: k=100, K1=2, K2=2
  // Time constant ≈ 1/(k*K1*K2/(1+K1+K1*K2)) = 7/(100*4) ≈ 0.018 s
  // A_gas0 = C/7 ≈ 1.43e-7
  auto r_loose = RunCascadeSystem(1e-3, 100.0, 2.0, 2.0);
  auto r_tight = RunCascadeSystem(1e-8, 100.0, 2.0, 2.0);

  // With the fix, tight balance atol should produce more steps
  EXPECT_GT(r_tight.accepted, r_loose.accepted) << "Tight atol should require more steps than loose atol. "
                                                << "loose=" << r_loose.accepted << " tight=" << r_tight.accepted;
}

/// @brief Verify that the algebraic balance variable does not go deeply negative.
///
/// Uses an ultra-stiff system (k=1e4, K1=100, K2=50) where the embedded error
/// estimate produces near-zero Yerror for all variables. Without the step-change
/// fix, the solver would accept a huge first step, overshooting P > C_total and
/// forcing A_gas deeply negative.
TEST(DAEAlgebraicError, AlgebraicVariableDoesNotOvershootDeeply)
{
  // Ultra-stiff: transient is ~5e-3 s, solver converges quickly
  auto r = RunCascadeSystem(1e-8, 1.0e4, 100.0, 50.0);

  // A_gas should stay non-negative (or very close to zero).
  // The analytical solution has A_gas >= 0 at all times.
  EXPECT_GE(r.min_A_gas, -1.0e-8) << "A_gas overshot deeply negative: min_A_gas=" << r.min_A_gas;

  // After 30s the system should be fully converted: P ≈ C_total, A_gas ≈ 0
  EXPECT_NEAR(r.P_final, 1.0e-6, 1.0e-8);
  EXPECT_NEAR(r.A_gas_final, 0.0, 1.0e-10);
}
