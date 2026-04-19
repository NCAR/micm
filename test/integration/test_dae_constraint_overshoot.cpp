// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression test for algebraic-variable overshoot in the Rosenbrock DAE solver.
//
// The NormalizedError() function previously excluded algebraic variables from
// the step-acceptance error norm.  This allowed the solver to accept large
// internal steps where differential species overshoot a conservation budget,
// forcing the algebraic "balance" variable negative — a physically impossible
// state that the continuous system can never reach.
//
// System:
//   Species: A (reactant), B (product), C (algebraic balance)
//   Kinetics: A -> B   (fast, rate k = 1e4)
//   Constraint: A + B + C = C_total   (C is algebraic)
//
// With enough initial A and a fast rate, the solver used to overshoot: B would
// exceed C_total, making C = C_total - A - B negative.  Including algebraic
// variables in the error norm causes the solver to reject those steps.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;

/// @brief Verify that conservation-constrained algebraic variables stay non-negative
///        when fast kinetics drain the pool.
TEST(DAEConstraintOvershoot, AlgebraicVariableStaysNonNegative)
{
  auto A = Species("A");
  auto B = Species("B");
  auto C = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

  // Fast reaction: A -> B with a large rate constant
  // This is analogous to SO2 oxidation producing SO4: fast enough that the
  // solver wants to take large steps.
  double k = 1.0e4;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                    .SetPhase(gas_phase)
                    .Build();

  // Conservation constraint: A + B + C = C_total
  // C is the algebraic variable (last in the terms list).
  // In the continuous system, C >= 0 always because A,B cannot exceed C_total
  // together. But the discrete solver can overshoot.
  double C_total = 1.0e-6;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint("mass_conservation", { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, C_total));

  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);

  std::size_t A_idx = state.variable_map_.at("A");
  std::size_t B_idx = state.variable_map_.at("B");
  std::size_t C_idx = state.variable_map_.at("C");

  // Initial conditions: most of the budget in A, a small amount in C, none in B.
  // The fast A->B reaction will rapidly convert A to B. The algebraic variable
  // C = C_total - A - B must remain >= 0 throughout.
  state.variables_[0][A_idx] = 0.9e-6;
  state.variables_[0][B_idx] = 0.0;
  state.variables_[0][C_idx] = 0.1e-6;  // C_total - A
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  // Integrate for 30 seconds with a large external time step.
  // The solver picks its own internal steps based on error control.
  double dt = 30.0;
  double advanced = 0.0;

  while (advanced < dt)
  {
    auto result = solver.Solve(dt - advanced, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Solver did not converge at t=" << advanced;
    advanced += result.stats_.final_time_;

    // Check conservation
    double sum = state.variables_[0][A_idx] + state.variables_[0][B_idx] + state.variables_[0][C_idx];
    EXPECT_NEAR(sum, C_total, 1.0e-12) << "Conservation violated at t=" << advanced;

    // The key assertion: C (the algebraic balance variable) must not go negative.
    // Before the fix, the solver would accept steps where B overshoots C_total,
    // causing C = C_total - A - B < 0.
    EXPECT_GE(state.variables_[0][C_idx], -1.0e-18)
        << "Algebraic variable C went negative (" << state.variables_[0][C_idx] << ") at t=" << advanced
        << "; A=" << state.variables_[0][A_idx] << ", B=" << state.variables_[0][B_idx];
  }

  // After 30s with k=1e4, A should be essentially 0 and B ~ A_init.
  // C retains its initial value since the reaction only converts A -> B.
  EXPECT_LT(state.variables_[0][A_idx], 1.0e-12);
  EXPECT_NEAR(state.variables_[0][B_idx], 0.9e-6, 1.0e-10);
  EXPECT_NEAR(state.variables_[0][C_idx], 0.1e-6, 1.0e-10);
}

/// @brief Same test with equilibrium + linear constraints (more species, closer to real cloud chemistry)
TEST(DAEConstraintOvershoot, EquilibriumPlusConservation)
{
  // System:
  //   A_gas  -- equilibrium -->  A_aq  (algebraic, K_eq * A_gas = A_aq)
  //   A_aq   -- fast kinetics -> P     (differential, rate = k * A_aq)
  //   Conservation: A_gas + A_aq + P = C_total   (A_gas is algebraic balance)
  //
  // The fast A_aq -> P reaction drains the pool. The equilibrium replenishes
  // A_aq from A_gas. The conservation constraint sets A_gas = C_total - A_aq - P.
  // If the solver overshoots P, A_gas goes negative.

  auto A_gas = Species("A_gas");
  auto A_aq = Species("A_aq");
  auto P = Species("P");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A_gas, A_aq, P } };

  // Fast reaction: A_aq -> P
  double k = 1.0e3;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A_aq })
                    .SetProducts({ { P, 1 } })
                    .SetRateConstant(ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                    .SetPhase(gas_phase)
                    .Build();

  double C_total = 1.0e-6;
  double K_eq = 10.0;

  std::vector<Constraint> constraints;

  // Equilibrium: K_eq * A_gas = A_aq  (A_aq is algebraic — first product)
  constraints.push_back(EquilibriumConstraint(
      "gas_aq_eq",
      std::vector<StoichSpecies>{ { A_gas, 1.0 } },
      std::vector<StoichSpecies>{ { A_aq, 1.0 } },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = 0.0 }));

  // Conservation: A_gas + A_aq + P = C_total  (A_gas is algebraic balance — last term)
  // Note: A_gas appears last so it becomes the algebraic variable for this constraint.
  constraints.push_back(LinearConstraint("mass_conservation", { { A_aq, 1.0 }, { P, 1.0 }, { A_gas, 1.0 } }, C_total));

  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);

  std::size_t A_gas_idx = state.variable_map_.at("A_gas");
  std::size_t A_aq_idx = state.variable_map_.at("A_aq");
  std::size_t P_idx = state.variable_map_.at("P");

  // Initial: most sulfur in gas phase, equilibrium satisfied, no product yet
  double A_gas_init = C_total / (1.0 + K_eq);  // ~ 9.09e-8
  double A_aq_init = K_eq * A_gas_init;         // ~ 9.09e-7
  state.variables_[0][A_gas_idx] = A_gas_init;
  state.variables_[0][A_aq_idx] = A_aq_init;
  state.variables_[0][P_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  double dt = 30.0;
  double advanced = 0.0;

  while (advanced < dt)
  {
    auto result = solver.Solve(dt - advanced, state);
    ASSERT_EQ(result.state_, SolverState::Converged) << "Solver did not converge at t=" << advanced;
    advanced += result.stats_.final_time_;

    double sum = state.variables_[0][A_gas_idx] + state.variables_[0][A_aq_idx] + state.variables_[0][P_idx];
    EXPECT_NEAR(sum, C_total, 1.0e-12) << "Conservation violated at t=" << advanced;

    // A_gas must not go negative
    EXPECT_GE(state.variables_[0][A_gas_idx], -1.0e-18)
        << "A_gas went negative (" << state.variables_[0][A_gas_idx] << ") at t=" << advanced;

    // A_aq must not go negative
    EXPECT_GE(state.variables_[0][A_aq_idx], -1.0e-18)
        << "A_aq went negative (" << state.variables_[0][A_aq_idx] << ") at t=" << advanced;
  }

  // After 30s with k=1e3, nearly all sulfur should be in P
  EXPECT_NEAR(state.variables_[0][P_idx], C_total, 1.0e-8);
  EXPECT_GE(state.variables_[0][A_gas_idx], -1.0e-18);
}
