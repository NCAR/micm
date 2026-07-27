// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <unordered_map>

using namespace micm;

/// @brief Helper: build a simple A→B system with equilibrium constraint C = K_eq * B
struct SimpleConstrainedSystem
{
  static constexpr double K = 0.5;
  static constexpr double K_EQ = 2.0;
  static constexpr double DELTA_H = 0.0;  // No temperature dependence for simplicity

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
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = K, .B_ = 0, .C_ = 0 })
                      .SetPhase(gas_phase)
                      .Build();

    // Equilibrium constraint: K_eq * B - C = 0, so C = K_eq * B
    std::vector<Constraint> constraints;
    constraints.emplace_back(EquilibriumConstraint(
        "B_C_eq",
        C,
        std::vector<StoichSpecies>{ { B, 1.0 } },
        std::vector<StoichSpecies>{ { C, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K_EQ, .delta_H_ = DELTA_H }));

    return builder.SetSystem(System(gas_phase))
        .SetReactions({ rxn })
        .SetConstraints(std::move(constraints))
        .SetReorderState(false)
        .Build();
  }
};

/// @brief Nonlinear constraint G = scale * (Z^2 - X) used to test row scaling and line search
class ScaledSquareRootConstraintModel
{
 public:
  explicit ScaledSquareRootConstraintModel(double row_scale)
      : row_scale_(row_scale)
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { "Z" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { "X", "Z" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& indices) const
  {
    const auto z = indices.at("Z");
    return { { z, indices.at("X") }, { z, z } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<Conditions>&, DenseMatrixPolicy&) { };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& indices) const
  {
    const auto x = indices.at("X");
    const auto z = indices.at("Z");
    const double row_scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& residual)
    {
      for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
      {
        residual[cell][z] = row_scale * (state[cell][z] * state[cell][z] - state[cell][x]);
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& indices,
      const SparseMatrixPolicy&) const
  {
    const auto x = indices.at("X");
    const auto z = indices.at("Z");
    const double row_scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
    {
      for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
      {
        jacobian[block][z][x] -= -row_scale;
        jacobian[block][z][z] -= 2.0 * row_scale * state[block][z];
      }
    };
  }

 private:
  double row_scale_;
};

/// @brief X'=-1 with algebraic copy Z=X, used to force post-solve clamping
class SlavedConstantForcingModel
{
 public:
  std::set<std::string> SpeciesUsed() const
  {
    return { "X", "Z" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<Conditions>&, DenseMatrixPolicy&) { };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& indices) const
  {
    const auto x = indices.at("X");
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
    {
      for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
      {
        forcing[cell][x] -= 1.0;
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>&,
      const SparseMatrixPolicy&) const
  {
    return [](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&) { };
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { "Z" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { "X", "Z" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& indices) const
  {
    const auto z = indices.at("Z");
    return { { z, indices.at("X") }, { z, z } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<Conditions>&, DenseMatrixPolicy&) { };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& indices) const
  {
    const auto x = indices.at("X");
    const auto z = indices.at("Z");
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& residual)
    {
      for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
      {
        residual[cell][z] = state[cell][z] - state[cell][x];
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& indices,
      const SparseMatrixPolicy&) const
  {
    const auto x = indices.at("X");
    const auto z = indices.at("Z");
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
    {
      for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
      {
        jacobian[block][z][x] -= -1.0;
        jacobian[block][z][z] -= 1.0;
      }
    };
  }
};

using StandardBuilder =
    CpuSolverBuilder<RosenbrockSolverParameters, Matrix<double>, SparseMatrix<double, SparseMatrixStandardOrdering>>;

auto BuildSquareRootSolver(const RosenbrockSolverParameters& parameters, double row_scale)
{
  const auto x = Species("X");
  const auto z = Species("Z");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ x, z } };
  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .AddExternalModel(ScaledSquareRootConstraintModel(row_scale))
      .SetReorderState(false)
      .Build();
}

auto BuildSlavedConstantForcingSolver(const RosenbrockSolverParameters& parameters)
{
  const auto x = Species("X");
  const auto z = Species("Z");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ x, z } };
  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .AddExternalModel(SlavedConstantForcingModel())
      .SetReorderState(false)
      .Build();
}

/// @brief Test that consistent initial conditions don't change state
TEST(ConstraintInitialization, ConsistentICsUnchanged)
{
  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_EQ * 0.5;  // C = K_eq * B = consistent
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;

  solver.UpdateStateParameters(state);

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
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
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
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
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
                    .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 0.5, .B_ = 0, .C_ = 0 })
                    .SetPhase(gas_phase)
                    .Build();

  auto options = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                    .SetSystem(System(gas_phase))
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
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(3);  // 3 grid cells

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  double K_eq = SimpleConstrainedSystem::K_EQ;

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
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(options));
  auto state = solver.GetState(1);

  auto A_idx = state.variable_map_.at("A");
  auto B_idx = state.variable_map_.at("B");
  auto C_idx = state.variable_map_.at("C");

  // Start consistent
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = SimpleConstrainedSystem::K_EQ * 0.5;
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

/// @brief A linear constraint solved by the last permitted update must report success
TEST(ConstraintInitialization, FinalAllowedNewtonUpdateReportsConvergence)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_iterations_ = 1;
  auto solver = SimpleConstrainedSystem::Build(StandardBuilder(parameters));
  auto state = solver.GetState(1);

  const auto A_idx = state.variable_map_.at("A");
  const auto B_idx = state.variable_map_.at("B");
  const auto C_idx = state.variable_map_.at("C");
  state.variables_[0][A_idx] = 1.0;
  state.variables_[0][B_idx] = 0.5;
  state.variables_[0][C_idx] = 5.0;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);

  EXPECT_EQ(status, SolverState::Converged);
  EXPECT_EQ(stats.constraint_init_iterations_, 1u);
  EXPECT_DOUBLE_EQ(state.variables_[0][A_idx], 1.0);
  EXPECT_DOUBLE_EQ(state.variables_[0][B_idx], 0.5);
  EXPECT_DOUBLE_EQ(state.variables_[0][C_idx], 1.0);
}

/// @brief Algebraically equivalent row scales must produce the same corrected state
TEST(ConstraintInitialization, ConvergenceIsInvariantToConstraintRowScaling)
{
  const std::array<double, 3> row_scales{ 1.0e-14, 1.0, 1.0e14 };
  for (const double row_scale : row_scales)
  {
    auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
    auto solver = BuildSquareRootSolver(parameters, row_scale);
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-8);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
    state.variables_[0][state.variable_map_.at("X")] = 1.0;
    state.variables_[0][state.variable_map_.at("Z")] = 2.0;

    SolverStats stats;
    const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);

    EXPECT_EQ(status, SolverState::Converged) << "row scale=" << row_scale;
    EXPECT_NEAR(state.variables_[0][state.variable_map_.at("Z")], 1.0, 1.0e-10) << "row scale=" << row_scale;
  }
}

/// @brief Backtracking globalizes Newton for a finite but difficult initial guess
TEST(ConstraintInitialization, BacktrackingConvergesNearZeroDerivative)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = BuildSquareRootSolver(parameters, 1.0);
  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-8);
  state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
  state.variables_[0][state.variable_map_.at("X")] = 1.0;
  state.variables_[0][state.variable_map_.at("Z")] = 1.0e-6;

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);

  EXPECT_EQ(status, SolverState::Converged);
  EXPECT_NEAR(state.variables_[0][state.variable_map_.at("Z")], 1.0, 1.0e-10);
  EXPECT_LE(stats.constraint_init_iterations_, parameters.constraint_init_max_iterations_);
}

/// @brief Failed initialization restores the complete caller-provided state
TEST(ConstraintInitialization, FailureDoesNotLeavePartialNewtonUpdate)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_iterations_ = 1;
  parameters.constraint_init_tolerance_ = 1.0e-12;
  auto solver = BuildSquareRootSolver(parameters, 1.0);
  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-8);
  state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
  const auto x = state.variable_map_.at("X");
  const auto z = state.variable_map_.at("Z");
  state.variables_[0][x] = 1.0;
  state.variables_[0][z] = 10.0;

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);

  EXPECT_EQ(status, SolverState::ConstraintInitializationFailed);
  EXPECT_DOUBLE_EQ(state.variables_[0][x], 1.0);
  EXPECT_DOUBLE_EQ(state.variables_[0][z], 10.0);
}

/// @brief Clamping a differential variable must be followed by algebraic projection
TEST(ConstraintInitialization, PostSolveClampReprojectsAlgebraicVariables)
{
  auto parameters = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = BuildSlavedConstantForcingSolver(parameters);
  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);
  state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-9));
  const auto x = state.variable_map_.at("X");
  const auto z = state.variable_map_.at("Z");
  state.variables_[0][x] = 0.5;
  state.variables_[0][z] = 0.5;

  const auto result = solver.Solve(1.0, state, parameters);

  EXPECT_EQ(result.state_, SolverState::Converged);
  EXPECT_DOUBLE_EQ(state.variables_[0][x], 0.0);
  EXPECT_DOUBLE_EQ(state.variables_[0][z], 0.0);
  EXPECT_EQ(result.stats_.constraint_init_iterations_, 2u);
}

/// @brief A usable partial solution must remain algebraically consistent after clamping
TEST(ConstraintInitialization, PartialSolveClampReprojectsWithoutChangingStatus)
{
  auto parameters = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  parameters.max_number_of_steps_ = 0;
  auto solver = BuildSlavedConstantForcingSolver(parameters);
  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);
  state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-9));
  const auto x = state.variable_map_.at("X");
  const auto z = state.variable_map_.at("Z");
  state.variables_[0][x] = 1.0e-9;
  state.variables_[0][z] = 1.0e-9;

  const auto result = solver.Solve(1.0, state, parameters);

  EXPECT_EQ(result.state_, SolverState::ConvergenceExceededMaxSteps);
  EXPECT_EQ(result.stats_.accepted_, 1u);
  EXPECT_DOUBLE_EQ(state.variables_[0][x], 0.0);
  EXPECT_DOUBLE_EQ(state.variables_[0][z], 0.0);
  EXPECT_EQ(result.stats_.constraint_init_iterations_, 2u);
}

/// @brief A wrapper-level initialization failure must not be altered by post-solve clamping
TEST(ConstraintInitialization, SolveFailurePreservesCallerState)
{
  auto parameters = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  parameters.constraint_init_max_iterations_ = 0;
  auto solver = BuildSlavedConstantForcingSolver(parameters);
  auto state = solver.GetState(1);
  const auto x = state.variable_map_.at("X");
  const auto z = state.variable_map_.at("Z");
  state.variables_[0][x] = -0.5;
  state.variables_[0][z] = -0.5;

  const auto result = solver.Solve(1.0, state, parameters);

  EXPECT_EQ(result.state_, SolverState::ConstraintInitializationFailed);
  EXPECT_DOUBLE_EQ(state.variables_[0][x], -0.5);
  EXPECT_DOUBLE_EQ(state.variables_[0][z], -0.5);
}
