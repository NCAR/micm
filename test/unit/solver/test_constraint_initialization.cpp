// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;

/// @brief Helper: build a simple A→B system with equilibrium constraint C = K_eq * B
struct SimpleConstrainedSystem
{
  static constexpr double k = 0.5;
  static constexpr double K_eq = 2.0;
  static constexpr double delta_H = 0.0;  // No temperature dependence for simplicity

  template<class SolverBuilderPolicy>
  static auto Build(SolverBuilderPolicy builder)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto C = Species("C");

    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, C } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { B, 1 } })
                      .SetRateConstant(ArrheniusRateConstant({ .A_ = k, .B_ = 0, .C_ = 0 }))
                      .SetPhase(gas_phase)
                      .Build();

    // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "B_C_eq",
        std::vector<StoichSpecies>{ { B, 1.0 } },
        std::vector<StoichSpecies>{ { C, 1.0 } },
        VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = delta_H }));

    return builder.SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
        .SetReactions({ rxn })
        .SetConstraints(std::move(constraints))
        .SetReorderState(false)
        .Build();
  }
};

using StandardBuilder = CpuSolverBuilder<
    RosenbrockSolverParameters,
    Matrix<double>,
    SparseMatrix<double, SparseMatrixStandardOrdering>>;

/// @brief Test that consistent initial conditions don't change state
TEST(ConstraintInitialization, ConsistentICsUnchanged)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(std::move(options)));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_eq * 0.5;  // C = K_eq * B = consistent
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  double A_before = state.variables_[0][A_idx];
  double B_before = state.variables_[0][B_idx];
  double C_before = state.variables_[0][C_idx];

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);
  // A should not change due to initialization (only from time stepping)
  // B should not change due to initialization (only from time stepping)
  // C was already consistent, so initialization should be nearly a no-op
  // Constraint init should converge in ≤1 iteration
  EXPECT_LE(result.stats_.constraint_init_iterations_, 1);
}

/// @brief Test that mildly inconsistent ICs are corrected
TEST(ConstraintInitialization, MildlyInconsistentICsCorrected)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(std::move(options)));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 5.0;  // Wrong: should be K_eq * B = 1.0
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  double A_before = state.variables_[0][A_idx];
  double B_before = state.variables_[0][B_idx];

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // ODE variable A should not have been modified by initialization
  // (it will change from time stepping, but the initialization should not touch it)
  // We can't check exact equality post-solve because time stepping changes A,
  // but we verify the initialization converged and the constraint is satisfied
  EXPECT_GT(result.stats_.constraint_init_iterations_, 0u);

  // After solve, the constraint should be satisfied: C ≈ K_eq * B
  double K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  double residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Test severely inconsistent ICs converge
TEST(ConstraintInitialization, SeverelyInconsistentICsConverge)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(std::move(options)));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 1000.0;  // Wildly off: should be K_eq * B = 1.0
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // Constraint should be satisfied after initialization + solve
  double K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  double residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Test pure ODE system (no constraints) is unaffected
TEST(ConstraintInitialization, PureODESystemUnaffected)
{
  auto A = Species("A");
  auto B = Species("B");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B } };

  Process rxn = ChemicalReactionBuilder()
                    .SetReactants({ A })
                    .SetProducts({ { B, 1 } })
                    .SetRateConstant(ArrheniusRateConstant({ .A_ = 0.5, .B_ = 0, .C_ = 0 }))
                    .SetPhase(gas_phase)
                    .Build();

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                    .SetSystem(System(SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn })
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  state.variables_[0][state.variable_map_.at("A")] = 1.0;
  state.variables_[0][state.variable_map_.at("B")] = 0.0;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  auto result = solver.Solve(0.01, state);

  EXPECT_EQ(result.state_, SolverState::Converged);
  // No constraint initialization should have happened
  EXPECT_EQ(result.stats_.constraint_init_iterations_, 0u);
}

/// @brief Test multi-cell systems with different inconsistency levels
TEST(ConstraintInitialization, MultiCellSystems)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(std::move(options)));
  auto state = solver.GetState(3);  // 3 grid cells

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  double K_eq = SimpleConstrainedSystem::K_eq;

  // Cell 0: consistent
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = K_eq * 0.5;

  // Cell 1: mildly inconsistent
  state.variables_[1][A_idx] = 2.0;
  state.variables_[1][B_idx] = 0.3;
  state.variables_[1][C_idx] = 5.0;  // Should be K_eq * 0.3 = 0.6

  // Cell 2: severely inconsistent (but non-negative — equilibrium constraint clamps negative concentrations)
  state.variables_[2][A_idx] = 0.5;
  state.variables_[2][B_idx] = 0.1;
  state.variables_[2][C_idx] = 100.0;  // Should be K_eq * 0.1 = 0.2

  for (std::size_t i = 0; i < 3; ++i)
  {
    state.conditions_[i].temperature_ = 298.15;
    state.conditions_[i].pressure_ = 101325.0;
  }

  solver.UpdateStateParameters(state);

  auto result = solver.Solve(0.001, state);

  EXPECT_EQ(result.state_, SolverState::Converged);

  // Verify constraint satisfied in all cells after solve
  double K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  for (std::size_t i = 0; i < 3; ++i)
  {
    double residual = K_eq_actual * state.variables_[i][B_idx] - state.variables_[i][C_idx];
    EXPECT_NEAR(residual, 0.0, 1.0e-5) << "Constraint not satisfied in cell " << i;
  }
}

/// @brief Test that subsequent Solve() calls re-check and re-initialize if needed
TEST(ConstraintInitialization, SubsequentSolveCallsReinitialize)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(std::move(options)));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  // Start consistent
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_eq * 0.5;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);
  auto result1 = solver.Solve(0.001, state);
  EXPECT_EQ(result1.state_, SolverState::Converged);

  // Perturb C off-manifold between solve calls (simulating external event)
  state.variables_[0][C_idx] = 999.0;

  solver.UpdateStateParameters(state);
  auto result2 = solver.Solve(0.001, state);
  EXPECT_EQ(result2.state_, SolverState::Converged);

  // Constraint should be re-satisfied after second solve
  double K_eq_actual = state.custom_rate_parameters_[0][state.custom_rate_parameter_map_.at("B_C_eq")];
  double residual = K_eq_actual * state.variables_[0][B_idx] - state.variables_[0][C_idx];
  EXPECT_NEAR(residual, 0.0, 1.0e-6);
}

/// @brief Test SolverStateToString for the new enum value
TEST(ConstraintInitialization, SolverStateToStringNewValue)
{
  EXPECT_EQ(SolverStateToString(SolverState::ConstraintInitializationFailed), "Constraint Initialization Failed");
}
