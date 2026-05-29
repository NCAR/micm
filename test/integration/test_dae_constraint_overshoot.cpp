// Copyright (c) 2023-2026 University Corporation for Atmospheric Research
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
//   Species: a (reactant), b (product), c (algebraic balance)
//   Kinetics: a -> b   (fast, rate k = 1e4)
//   Constraint: a + b + c = c_total   (c is algebraic)
//
// With enough initial a and a fast rate, the solver used to overshoot: b would
// exceed c_total, making c = c_total - a - b negative.  Including algebraic
// variables in the error norm causes the solver to reject those steps.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <string>

using namespace micm;

/// @brief Verify that conservation-constrained algebraic variables stay non-negative
///        when fast kinetics drain the pool.
TEST(DAEConstraintOvershoot, AlgebraicVariableStaysNonNegative)
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a, b, c } };

  // Fast reaction: a -> b with a large rate constant
  // This is analogous to SO2 oxidation producing SO4: fast enough that the
  // solver wants to take large steps.
  double k = 1.0e4;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ a })
                    .SetProducts({ { b, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  // Conservation constraint: a + b + c = c_total
  // c is the explicitly set algebraic variable.
  // In the continuous system, c >= 0 always because a,b cannot exceed c_total
  // together. But the discrete solver can overshoot.
  double c_total = 1.0e-6;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint("mass_conservation", c, { { a, 1.0 }, { b, 1.0 }, { c, 1.0 } }, c_total));

  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);
  state.SetAbsoluteTolerances(std::vector<double>(3, 1.0e-12));

  std::size_t a_idx = state.variable_map_.at("A");
  std::size_t b_idx = state.variable_map_.at("B");
  std::size_t c_idx = state.variable_map_.at("C");

  // Initial conditions: most of the budget in a, a small amount in c, none in b.
  // The fast a->b reaction will rapidly convert a to b. The algebraic variable
  // c = c_total - a - b must remain >= 0 throughout.
  state.variables_[0][a_idx] = 0.9e-6;
  state.variables_[0][b_idx] = 0.0;
  state.variables_[0][c_idx] = 0.1e-6;  // c_total - a
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
    ASSERT_EQ(result.state_, SolverState::CONVERGED) << "Solver did not converge at t=" << advanced;
    advanced += result.stats_.final_time_;

    // Check conservation
    double sum = state.variables_[0][a_idx] + state.variables_[0][b_idx] + state.variables_[0][c_idx];
    EXPECT_NEAR(sum, c_total, 1.0e-12) << "Conservation violated at t=" << advanced;

    // The key assertion: c (the algebraic balance variable) must not go negative.
    // Before the fix, the solver would accept steps where b overshoots c_total,
    // causing c = c_total - a - b < 0.
    EXPECT_GE(state.variables_[0][c_idx], -1.0e-18)
        << "Algebraic variable C went negative (" << state.variables_[0][c_idx] << ") at t=" << advanced
        << "; A=" << state.variables_[0][a_idx] << ", B=" << state.variables_[0][b_idx];
  }

  // After 30s with k=1e4, a should be essentially 0 and b ~ A_init.
  // c retains its initial value since the reaction only converts a -> b.
  EXPECT_LT(state.variables_[0][a_idx], 1.0e-12);
  EXPECT_NEAR(state.variables_[0][b_idx], 0.9e-6, 1.0e-10);
  EXPECT_NEAR(state.variables_[0][c_idx], 0.1e-6, 1.0e-10);
}

/// @brief Same test with equilibrium + linear constraints (more species, closer to real cloud chemistry)
TEST(DAEConstraintOvershoot, EquilibriumPlusConservation)
{
  // System:
  //   a_gas  -- equilibrium -->  a_aq  (algebraic, K_eq * a_gas = a_aq)
  //   a_aq   -- fast kinetics -> p     (differential, rate = k * a_aq)
  //   Conservation: a_gas + a_aq + p = c_total   (a_gas is algebraic balance)
  //
  // The fast a_aq -> p reaction drains the pool. The equilibrium replenishes
  // a_aq from a_gas. The conservation constraint sets a_gas = c_total - a_aq - p.
  // If the solver overshoots p, a_gas goes negative.

  auto a_gas = Species("A_gas");
  auto a_aq = Species("A_aq");
  auto p = Species("P");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a_gas, a_aq, p } };

  // Fast reaction: a_aq -> p
  double k = 1.0e3;
  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ a_aq })
                    .SetProducts({ { p, 1 } })
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  double c_total = 1.0e-6;
  double k_eq = 10.0;

  std::vector<Constraint> constraints;

  // Equilibrium: K_eq * a_gas = a_aq  (a_aq is the explicitly set algebraic species)
  constraints.push_back(EquilibriumConstraint(
      "gas_aq_eq",
      a_aq,
      std::vector<StoichSpecies>{ { a_gas, 1.0 } },
      std::vector<StoichSpecies>{ { a_aq, 1.0 } },
      VantHoffParam{ .K_HLC_ref_ = k_eq, .delta_H_ = 0.0 }));

  // Conservation: a_gas + a_aq + p = c_total  (a_gas is the algebraic balance variable)
  constraints.push_back(
      LinearConstraint("mass_conservation", a_gas, { { a_aq, 1.0 }, { p, 1.0 }, { a_gas, 1.0 } }, c_total));

  auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(gas_phase))
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);

  std::size_t a_gas_idx = state.variable_map_.at("A_gas");
  std::size_t a_aq_idx = state.variable_map_.at("A_aq");
  std::size_t p_idx = state.variable_map_.at("P");

  // Use reasonable absolute tolerances:
  // - Differential variable (p): tight tolerance for accuracy
  // - Algebraic variables (a_gas, a_aq): moderate tolerance to allow legitimate step changes
  //   while still detecting overshoot via the step-change error estimate
  std::vector<double> atols(3, 1.0e-12);
  atols[a_gas_idx] = 1.0e-8;
  atols[a_aq_idx] = 1.0e-8;
  state.SetAbsoluteTolerances(atols);

  // Initial: most sulfur in gas phase, equilibrium satisfied, no product yet
  double a_gas_init = c_total / (1.0 + k_eq);  // ~ 9.09e-8
  double a_aq_init = k_eq * a_gas_init;        // ~ 9.09e-7
  state.variables_[0][a_gas_idx] = a_gas_init;
  state.variables_[0][a_aq_idx] = a_aq_init;
  state.variables_[0][p_idx] = 0.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  double dt = 30.0;
  double advanced = 0.0;

  while (advanced < dt)
  {
    auto result = solver.Solve(dt - advanced, state);
    ASSERT_EQ(result.state_, SolverState::CONVERGED) << "Solver did not converge at t=" << advanced;
    advanced += result.stats_.final_time_;

    double sum = state.variables_[0][a_gas_idx] + state.variables_[0][a_aq_idx] + state.variables_[0][p_idx];
    EXPECT_NEAR(sum, c_total, 1.0e-12) << "Conservation violated at t=" << advanced;

    // a_gas must not go negative
    EXPECT_GE(state.variables_[0][a_gas_idx], -1.0e-18)
        << "A_gas went negative (" << state.variables_[0][a_gas_idx] << ") at t=" << advanced;

    // a_aq must not go negative
    EXPECT_GE(state.variables_[0][a_aq_idx], -1.0e-18)
        << "A_aq went negative (" << state.variables_[0][a_aq_idx] << ") at t=" << advanced;
  }

  // After 30s with k=1e3, nearly all sulfur should be in p
  EXPECT_NEAR(state.variables_[0][p_idx], c_total, 1.0e-8);
  EXPECT_GE(state.variables_[0][a_gas_idx], -1.0e-18);
}

/// @brief Exercise conservation-constrained overshoot behavior for all Rosenbrock parameter sets.
TEST(DAEConstraintOvershoot, AllRosenbrockOrdersConstrained)
{
  // The two-stage Rosenbrock lacks the stability to handle highly-stiff constrained systems
  // (k=1e4), so it is omitted here. Its normalized-error behavior with constraints is covered
  // by RosenbrockSolver.StandardNormalizedErrorWithConstraints.
  const std::vector<std::pair<std::string, RosenbrockSolverParameters>> parameter_sets = {
    { "three-stage", RosenbrockSolverParameters::ThreeStageRosenbrockParameters() },
    { "four-stage", RosenbrockSolverParameters::FourStageRosenbrockParameters() },
    { "four-stage-dae", RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters() },
    { "six-stage-dae", RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters() },
  };

  for (const auto& [name, options] : parameter_sets)
  {
    SCOPED_TRACE(name);

    auto a = Species("A");
    auto b = Species("B");
    auto c = Species("C");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a, b, c } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ { b, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0e4, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    constexpr double c_total = 1.0e-6;
    std::vector<Constraint> constraints;
    constraints.push_back(LinearConstraint("mass_conservation", c, { { a, 1.0 }, { b, 1.0 }, { c, 1.0 } }, c_total));

    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(3, 1.0e-12));

    const std::size_t a_idx = state.variable_map_.at("A");
    const std::size_t b_idx = state.variable_map_.at("B");
    const std::size_t c_idx = state.variable_map_.at("C");

    state.variables_[0][a_idx] = 0.9e-6;
    state.variables_[0][b_idx] = 0.0;
    state.variables_[0][c_idx] = 0.1e-6;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;

    solver.UpdateStateParameters(state);

    double advanced = 0.0;
    constexpr double dt = 30.0;
    while (advanced < dt)
    {
      auto result = solver.Solve(dt - advanced, state);
      ASSERT_EQ(result.state_, SolverState::CONVERGED) << "Solver did not converge for " << name << " at t=" << advanced;
      advanced += result.stats_.final_time_;

      const double sum = state.variables_[0][a_idx] + state.variables_[0][b_idx] + state.variables_[0][c_idx];
      EXPECT_NEAR(sum, c_total, 1.0e-12);
      EXPECT_GE(state.variables_[0][c_idx], -1.0e-18)
          << "Algebraic variable C went negative for " << name << ": " << state.variables_[0][c_idx];
    }

    EXPECT_LT(state.variables_[0][a_idx], 1.0e-12);
    EXPECT_NEAR(state.variables_[0][b_idx], 0.9e-6, 1.0e-10);
    EXPECT_NEAR(state.variables_[0][c_idx], 0.1e-6, 1.0e-10);
  }
}
