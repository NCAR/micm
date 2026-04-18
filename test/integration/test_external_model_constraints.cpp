// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "stub_aerosol_with_constraints.hpp"

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/util/jacobian_verification.hpp>

#include <gtest/gtest.h>

#include <cmath>

/// @brief Constraint-only external model that enforces K_eq * [reactant] - [product] = 0
///
/// This model contributes no state variables and no processes — only an algebraic
/// equilibrium constraint. Used with AddExternalModel() to test that
/// external model constraints produce results equivalent to kinetic systems and
/// built-in EquilibriumConstraint.
class EquilibriumConstraintModel
{
 public:
  EquilibriumConstraintModel(const std::string& reactant, const std::string& product, double K_eq)
      : reactant_(reactant),
        product_(product),
        K_eq_(K_eq)
  {
  }

  // ──── Constraint methods (satisfies HasConstraints concept) ────

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { product_ };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { reactant_, product_ };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    auto i_r = state_indices.at(reactant_);
    auto i_p = state_indices.at(product_);
    return { { i_p, i_r }, { i_p, i_p } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  /// Residual: G = K_eq * [reactant] - [product]
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var) const
  {
    auto i_r = var.at(reactant_);
    auto i_p = var.at(product_);
    double K = K_eq_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
        forcing[i][i_p] = K * state[i][i_r] - state[i][i_p];
    };
  }

  /// Jacobian: dG/d[reactant] = K_eq, dG/d[product] = -1
  /// Subtracted per solver convention: jac -= dG/dy
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var,
      const SparseMatrixPolicy&) const
  {
    auto i_r = var.at(reactant_);
    auto i_p = var.at(product_);
    double K = K_eq_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
      {
        jac[i][i_p][i_r] -= K;     // -dG/d[reactant] = -K_eq
        jac[i][i_p][i_p] -= -1.0;  // -dG/d[product]  = +1
      }
    };
  }

 private:
  std::string reactant_;
  std::string product_;
  double K_eq_;
};

/// @brief Constraint-only external model providing BOTH equilibrium AND conservation constraints
///
/// For a system A → B with fast B ⇌ C equilibrium:
///   - B row (conservation):  G_B = [A] + [B] + [C] - total = 0
///   - C row (equilibrium):   G_C = K_eq * [B] - [C] = 0
///
/// Both B and C are algebraic. Only A evolves via ODE.
/// This makes the constraint system equivalent to the kinetic system in steady state:
///   [A]=0, [B]=total/(1+K_eq), [C]=K_eq*total/(1+K_eq)
class ConservativeEquilibriumConstraintModel
{
 public:
  ConservativeEquilibriumConstraintModel(
      const std::string& a,
      const std::string& b,
      const std::string& c,
      double K_eq,
      double total)
      : species_a_(a),
        species_b_(b),
        species_c_(c),
        K_eq_(K_eq),
        total_(total)
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { species_b_, species_c_ };  // both algebraic
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { species_a_, species_b_, species_c_ };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    auto i_a = state_indices.at(species_a_);
    auto i_b = state_indices.at(species_b_);
    auto i_c = state_indices.at(species_c_);
    return {
      // B row (conservation): dG_B/dA=1, dG_B/dB=1, dG_B/dC=1
      { i_b, i_a }, { i_b, i_b }, { i_b, i_c },
      // C row (equilibrium): dG_C/dB=K_eq, dG_C/dC=-1
      { i_c, i_b }, { i_c, i_c }
    };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var) const
  {
    auto i_a = var.at(species_a_);
    auto i_b = var.at(species_b_);
    auto i_c = var.at(species_c_);
    double K = K_eq_;
    double tot = total_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        // B row: conservation
        forcing[i][i_b] = state[i][i_a] + state[i][i_b] + state[i][i_c] - tot;
        // C row: equilibrium
        forcing[i][i_c] = K * state[i][i_b] - state[i][i_c];
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var,
      const SparseMatrixPolicy&) const
  {
    auto i_a = var.at(species_a_);
    auto i_b = var.at(species_b_);
    auto i_c = var.at(species_c_);
    double K = K_eq_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
      {
        // B row (conservation): dG_B/dA=1, dG_B/dB=1, dG_B/dC=1
        jac[i][i_b][i_a] -= 1.0;
        jac[i][i_b][i_b] -= 1.0;
        jac[i][i_b][i_c] -= 1.0;
        // C row (equilibrium): dG_C/dB=K_eq, dG_C/dC=-1
        jac[i][i_c][i_b] -= K;
        jac[i][i_c][i_c] -= -1.0;
      }
    };
  }

 private:
  std::string species_a_;
  std::string species_b_;
  std::string species_c_;
  double K_eq_;
  double total_;
};

/// @brief Constraint-only external model that enforces mass conservation for one species pool
///
/// Replaces the ODE row of `controlled_species_` with:
///   G = sum([all species]) - total = 0
///
/// This is used alongside EquilibriumConstraintModels to compose conservative equilibrium
/// systems from separate models (equilibrium ratio + mass conservation).
class MassConservationModel
{
 public:
  MassConservationModel(const std::string& controlled_species, const std::vector<std::string>& all_species, double total)
      : controlled_species_(controlled_species),
        all_species_(all_species),
        total_(total)
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { controlled_species_ };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { all_species_.begin(), all_species_.end() };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    auto i_ctrl = state_indices.at(controlled_species_);
    std::set<std::pair<std::size_t, std::size_t>> elements;
    for (auto& sp : all_species_)
      elements.insert({ i_ctrl, state_indices.at(sp) });
    return elements;
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var) const
  {
    auto i_ctrl = var.at(controlled_species_);
    std::vector<std::size_t> indices;
    for (auto& sp : all_species_)
      indices.push_back(var.at(sp));
    double tot = total_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        double sum = 0.0;
        for (auto idx : indices)
          sum += state[i][idx];
        forcing[i][i_ctrl] = sum - tot;
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& var,
      const SparseMatrixPolicy&) const
  {
    auto i_ctrl = var.at(controlled_species_);
    std::vector<std::size_t> indices;
    for (auto& sp : all_species_)
      indices.push_back(var.at(sp));
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
        for (auto idx : indices)
          jac[i][i_ctrl][idx] -= 1.0;  // dG/d[species] = 1
    };
  }

 private:
  std::string controlled_species_;
  std::vector<std::string> all_species_;
  double total_;
};

/// @brief Verify that AddExternalModel wraps both processes and constraints for a constrained model
TEST(ExternalModelConstraints, AddExternalModelWithConstraints)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify constraint metadata
  EXPECT_EQ(state.state_size_, 2);       // A_GAS, AEROSOL.A_AQ
  EXPECT_EQ(state.constraint_size_, 1);  // one algebraic constraint

  // Verify mass matrix diagonal: A_AQ row should be algebraic (0.0)
  auto i_gas = state.variable_map_.at("A_GAS");
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_gas], 1.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);
}

/// @brief Verify that AddExternalModel works for a process-only model (no constraints)
TEST(ExternalModelConstraints, AddExternalModelProcessOnly)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  // No total_mass → constraints disabled
  StubAerosolWithConstraints aerosol(0.01);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // No constraints — pure ODE mode
  EXPECT_EQ(state.constraint_size_, 0);
  auto i_gas = state.variable_map_.at("A_GAS");
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_gas], 1.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 1.0);
}

/// @brief Verify that the DAE solver enforces the mass conservation constraint
TEST(ExternalModelConstraints, DAESolveEnforcesConservation)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  double k = 0.1;
  StubAerosolWithConstraints aerosol(k, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Initialize: most mass in gas phase
  state.variables_[0][state.variable_map_.at("A_GAS")] = 0.9;
  state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")] = 0.1;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Solve several time steps and verify conservation
  double dt = 10.0;
  for (int step = 0; step < 20; ++step)
  {
    auto result = solver.Solve(dt, state);
    EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Step " << step;

    double sum = state.variables_[0][state.variable_map_.at("A_GAS")]
               + state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")];
    EXPECT_NEAR(sum, total, 1e-4) << "Conservation violated at step " << step;
  }
}

/// @brief Verify that external model constraints combine with built-in SetConstraints
TEST(ExternalModelConstraints, CombinedBuiltInAndExternalConstraints)
{
  auto A_GAS = micm::Species("A_GAS");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  micm::Phase gas_phase{ "gas", { A_GAS, B, C } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  // Built-in constraint: B <-> C equilibrium
  double K_eq = 5.0;
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq",
      std::vector<micm::StoichSpecies>{ { B, 1.0 } },
      std::vector<micm::StoichSpecies>{ { C, 1.0 } },
      micm::VantHoffParam{ K_eq, 0.0 }));

  // Process: A_GAS -> B
  double k_rxn = 0.05;
  micm::Process rxn = micm::ChemicalReactionBuilder()
                          .SetReactants({ A_GAS })
                          .SetProducts({ { B, 1 } })
                          .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_rxn, .B_ = 0, .C_ = 0 }))
                          .SetPhase(gas_phase)
                          .Build();

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({ rxn })
                    .SetConstraints(std::move(constraints))
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // 2 algebraic variables: C (from built-in) and AEROSOL.A_AQ (from external)
  EXPECT_EQ(state.constraint_size_, 2);

  // Verify mass matrix diagonal
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  auto i_c = state.variable_map_.at("C");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_c], 0.0);
}

/// @brief Verify AddExternalModel() adds only processes when using standard Rosenbrock parameters
TEST(ExternalModelConstraints, AddExternalModelOnlyStandardRosenbrock)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  StubAerosolWithConstraints aerosol(0.01);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // Add only processes (constraints are not enabled with standard Rosenbrock parameters)
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);
  EXPECT_EQ(state.constraint_size_, 0);
}

/// @brief Verify AddExternalModel() adds only constraints when model lacks processes
TEST(ExternalModelConstraints, AddExternalModelConstraintsOnly)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.01, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Constraint should be active
  EXPECT_EQ(state.constraint_size_, 1);
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);

  // Solve and verify conservation
  state.variables_[0][state.variable_map_.at("A_GAS")] = 0.8;
  state.variables_[0][i_aq] = 0.2;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  auto result = solver.Solve(50.0, state);
  EXPECT_EQ(result.state_, micm::SolverState::Converged);

  double sum = state.variables_[0][state.variable_map_.at("A_GAS")]
             + state.variables_[0][i_aq];
  EXPECT_NEAR(sum, total, 1e-4);
}

/// @brief Verify multiple grid cells work with external model constraints
TEST(ExternalModelConstraints, MultiGridCell)
{
  auto A_GAS = micm::Species("A_GAS");
  micm::Phase gas_phase{ "gas", { A_GAS } };

  double total = 1.0;
  StubAerosolWithConstraints aerosol(0.1, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  const int num_cells = 3;
  auto state = solver.GetState(num_cells);

  // Different initial conditions per cell
  for (int c = 0; c < num_cells; ++c)
  {
    double gas_frac = 0.9 - 0.2 * c;
    state.variables_[c][state.variable_map_.at("A_GAS")] = gas_frac;
    state.variables_[c][state.variable_map_.at("AEROSOL.A_AQ")] = total - gas_frac;
    state.conditions_[c].temperature_ = 298.0;
    state.conditions_[c].pressure_ = 101325.0;
  }

  auto result = solver.Solve(50.0, state);
  EXPECT_EQ(result.state_, micm::SolverState::Converged);

  for (int c = 0; c < num_cells; ++c)
  {
    double sum = state.variables_[c][state.variable_map_.at("A_GAS")]
               + state.variables_[c][state.variable_map_.at("AEROSOL.A_AQ")];
    EXPECT_NEAR(sum, total, 1e-4) << "Conservation violated in cell " << c;
  }
}

// ───────────────────────────────────────────────────────────────────────────────
// Convergence tests: kinetic (process-only) vs constraint-based systems
//
// A single equilibrium constraint K_eq*[B]-[C]=0 enforces the equilibrium ratio
// but does NOT conserve mass (C is created algebraically from B's value).
// To get the same absolute concentrations as the kinetic system, the constraint
// formulation must include BOTH equilibrium AND mass conservation constraints,
// making B and C algebraic while A remains the only ODE variable.
//
// System: A → B (slow driver, k=0.1), B ⇌ C (fast, K_eq=5)
// Kinetic (ODE):       A→B (k), B→C (k_f), C→B (k_b)  — mass conserved by construction
// Constraint (DAE):    A→B (k), conservation + equilibrium constraints
// Steady state: [A]≈0, [B]=1/(1+K_eq), [C]=K_eq/(1+K_eq)
// ───────────────────────────────────────────────────────────────────────────────

namespace
{
  // Shared parameters for convergence tests
  constexpr double K_EQ = 5.0;       // equilibrium constant [C]/[B]
  constexpr double K_DRIVE = 0.1;    // rate constant for A → B (slow driver)
  constexpr double K_FWD = 10.0;     // rate constant for B → C (fast)
  constexpr double K_BWD = K_FWD / K_EQ;  // rate constant for C → B

  // Expected steady state for total=1.0
  constexpr double EXPECTED_A = 0.0;
  constexpr double EXPECTED_B = 1.0 / (1.0 + K_EQ);
  constexpr double EXPECTED_C = K_EQ / (1.0 + K_EQ);

  /// Helper: build and solve a kinetic (ODE) system with forward+backward reactions
  /// Returns (final_A, final_B, final_C)
  std::tuple<double, double, double> SolveKineticSystem()
  {
    auto A = micm::Species("A");
    auto B = micm::Species("B");
    auto C = micm::Species("C");
    micm::Phase gas_phase{ "gas", { A, B, C } };

    micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                               .SetReactants({ A })
                               .SetProducts({ { B, 1 } })
                               .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_DRIVE, .B_ = 0, .C_ = 0 }))
                               .SetPhase(gas_phase)
                               .Build();
    micm::Process rxn_bc = micm::ChemicalReactionBuilder()
                               .SetReactants({ B })
                               .SetProducts({ { C, 1 } })
                               .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_FWD, .B_ = 0, .C_ = 0 }))
                               .SetPhase(gas_phase)
                               .Build();
    micm::Process rxn_cb = micm::ChemicalReactionBuilder()
                               .SetReactants({ C })
                               .SetProducts({ { B, 1 } })
                               .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_BWD, .B_ = 0, .C_ = 0 }))
                               .SetPhase(gas_phase)
                               .Build();

    auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                      .SetReactions({ rxn_ab, rxn_bc, rxn_cb })
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.variables_[0][state.variable_map_.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;

    // Integrate well past all timescales:
    //   A→B timescale ~1/k = 10, B⇌C timescale ~1/(k_f+k_b) ≈ 0.08
    //   Total integration: 200 (20× the slowest timescale)
    double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Kinetic solve failed at step " << step;
    }

    return { state.variables_[0][state.variable_map_.at("A")],
             state.variables_[0][state.variable_map_.at("B")],
             state.variables_[0][state.variable_map_.at("C")] };
  }

  /// Helper: build and solve a DAE system using ConservativeEquilibriumConstraintModel
  /// (equilibrium + conservation constraints, making B and C algebraic)
  /// Returns (final_A, final_B, final_C)
  std::tuple<double, double, double> SolveConservativeConstraintSystem()
  {
    auto A = micm::Species("A");
    auto B = micm::Species("B");
    auto C = micm::Species("C");
    micm::Phase gas_phase{ "gas", { A, B, C } };

    // Only the slow driver — equilibrium is captured by constraints
    micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                               .SetReactants({ A })
                               .SetProducts({ { B, 1 } })
                               .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_DRIVE, .B_ = 0, .C_ = 0 }))
                               .SetPhase(gas_phase)
                               .Build();

    // Conservation + equilibrium constraints (B and C are algebraic)
    ConservativeEquilibriumConstraintModel eq_model("A", "B", "C", K_EQ, 1.0);

    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                      .SetReactions({ rxn_ab })
                      .AddExternalModel(eq_model)
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.variables_[0][state.variable_map_.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;

    double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Conservative constraint solve failed at step " << step;
    }

    return { state.variables_[0][state.variable_map_.at("A")],
             state.variables_[0][state.variable_map_.at("B")],
             state.variables_[0][state.variable_map_.at("C")] };
  }

  /// Helper: build and solve a DAE system using the simple EquilibriumConstraintModel
  /// (ratio constraint only — does NOT conserve mass)
  /// Returns (final_A, final_B, final_C)
  std::tuple<double, double, double> SolveSimpleConstraintSystem()
  {
    auto A = micm::Species("A");
    auto B = micm::Species("B");
    auto C = micm::Species("C");
    micm::Phase gas_phase{ "gas", { A, B, C } };

    micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                               .SetReactants({ A })
                               .SetProducts({ { B, 1 } })
                               .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_DRIVE, .B_ = 0, .C_ = 0 }))
                               .SetPhase(gas_phase)
                               .Build();

    // Simple ratio constraint: K_eq*[B] - [C] = 0 (C is algebraic, no conservation)
    EquilibriumConstraintModel eq_model("B", "C", K_EQ);

    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                      .SetReactions({ rxn_ab })
                      .AddExternalModel(eq_model)
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.variables_[0][state.variable_map_.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;

    double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Simple constraint solve failed at step " << step;
    }

    return { state.variables_[0][state.variable_map_.at("A")],
             state.variables_[0][state.variable_map_.at("B")],
             state.variables_[0][state.variable_map_.at("C")] };
  }
}  // anonymous namespace

/// @brief Kinetic and conservative constraint systems converge to the same steady state.
///
/// The conservative constraint system uses both equilibrium (K_eq*[B]-[C]=0) and mass
/// conservation ([A]+[B]+[C]-total=0) constraints, making B and C algebraic. This is
/// the correct DAE formulation equivalent to the kinetic B⇌C exchange: both distribute
/// mass identically at steady state.
///
/// System: A → B (k=0.1), B ⇌ C (K_eq=5)
/// Integration: 200s at dt=1s (20× the slowest timescale 1/0.1=10s)
/// Expected: [A]≈0, [B]≈1/6≈0.1667, [C]≈5/6≈0.8333
TEST(ExternalModelConstraints, KineticVsConservativeConstraintConvergence)
{
  auto [kin_A, kin_B, kin_C] = SolveKineticSystem();
  auto [con_A, con_B, con_C] = SolveConservativeConstraintSystem();

  // Both near analytical steady state
  EXPECT_NEAR(kin_A, EXPECTED_A, 1e-6);
  EXPECT_NEAR(con_A, EXPECTED_A, 1e-6);
  EXPECT_NEAR(kin_B, EXPECTED_B, 1e-3);
  EXPECT_NEAR(con_B, EXPECTED_B, 1e-3);
  EXPECT_NEAR(kin_C, EXPECTED_C, 1e-3);
  EXPECT_NEAR(con_C, EXPECTED_C, 1e-3);

  // The two systems agree with each other
  EXPECT_NEAR(kin_A, con_A, 1e-6);
  EXPECT_NEAR(kin_B, con_B, 1e-3);
  EXPECT_NEAR(kin_C, con_C, 1e-3);

  // Mass conservation in both
  EXPECT_NEAR(kin_A + kin_B + kin_C, 1.0, 1e-6);
  EXPECT_NEAR(con_A + con_B + con_C, 1.0, 1e-3);
}

/// @brief Simple equilibrium constraint preserves the ratio C/B = K_eq, but does NOT
///        conserve total mass. The kinetic and constraint systems should agree on the
///        equilibrium ratio, even though their absolute concentrations differ.
///
/// Without a conservation constraint, the algebraic C=K_eq*B creates mass: at steady
/// state [A]+[B]+[C] = 1+K_eq rather than 1.0. This is the expected behaviour for a
/// pure ratio constraint.
TEST(ExternalModelConstraints, SimpleConstraintPreservesRatioNotMass)
{
  auto [kin_A, kin_B, kin_C] = SolveKineticSystem();
  auto [sim_A, sim_B, sim_C] = SolveSimpleConstraintSystem();

  // Both enforce the equilibrium ratio
  EXPECT_NEAR(kin_C / kin_B, K_EQ, 1e-3);
  EXPECT_NEAR(sim_C / sim_B, K_EQ, 1e-3);

  // Kinetic system conserves mass
  EXPECT_NEAR(kin_A + kin_B + kin_C, 1.0, 1e-6);

  // Simple constraint system does NOT conserve total mass — C is created algebraically.
  // At t→∞: A≈0, B≈1 (all A→B mass), C≈K_eq*1=5, so total≈1+K_eq=6
  EXPECT_NEAR(sim_A + sim_B + sim_C, 1.0 + K_EQ, 0.1);
}

/// @brief Built-in EquilibriumConstraint and external model constraint are both DAE
///        systems and should produce nearly identical results at every time step.
TEST(ExternalModelConstraints, BuiltInVsExternalModelConstraintStepByStep)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  micm::Phase gas_phase{ "gas", { A, B, C } };

  micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                             .SetReactants({ A })
                             .SetProducts({ { B, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_DRIVE, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();

  // Built-in constraint solver
  std::vector<micm::Constraint> constraints;
  constraints.push_back(micm::EquilibriumConstraint(
      "B_C_eq", std::vector<micm::StoichSpecies>{ { B, 1.0 } }, std::vector<micm::StoichSpecies>{ { C, 1.0 } }, micm::VantHoffParam{ K_EQ, 0.0 }));

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto builtin_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                            .SetReactions({ rxn_ab })
                            .SetConstraints(std::move(constraints))
                            .SetReorderState(false)
                            .Build();

  // External model constraint solver
  EquilibriumConstraintModel eq_model("B", "C", K_EQ);
  auto ext_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                        .SetReactions({ rxn_ab })
                        .AddExternalModel(eq_model)
                        .SetReorderState(false)
                        .Build();

  auto state_bi = builtin_solver.GetState(1);
  auto state_ext = ext_solver.GetState(1);

  for (auto* s : { &state_bi, &state_ext })
  {
    s->variables_[0][s->variable_map_.at("A")] = 1.0;
    s->variables_[0][s->variable_map_.at("B")] = 0.0;
    s->variables_[0][s->variable_map_.at("C")] = 0.0;
    s->conditions_[0].temperature_ = 298.0;
    s->conditions_[0].pressure_ = 101325.0;
  }

  double dt = 0.5;
  for (int step = 0; step < 100; ++step)
  {
    builtin_solver.UpdateStateParameters(state_bi);
    ext_solver.UpdateStateParameters(state_ext);

    auto res_bi = builtin_solver.Solve(dt, state_bi);
    auto res_ext = ext_solver.Solve(dt, state_ext);

    ASSERT_EQ(res_bi.state_, micm::SolverState::Converged) << "Built-in failed step " << step;
    ASSERT_EQ(res_ext.state_, micm::SolverState::Converged) << "External failed step " << step;

    double bi_A = state_bi.variables_[0][state_bi.variable_map_.at("A")];
    double bi_B = state_bi.variables_[0][state_bi.variable_map_.at("B")];
    double bi_C = state_bi.variables_[0][state_bi.variable_map_.at("C")];
    double ext_A_val = state_ext.variables_[0][state_ext.variable_map_.at("A")];
    double ext_B_val = state_ext.variables_[0][state_ext.variable_map_.at("B")];
    double ext_C_val = state_ext.variables_[0][state_ext.variable_map_.at("C")];

    EXPECT_NEAR(bi_A, ext_A_val, 1e-8) << "A diverged at step " << step;
    EXPECT_NEAR(bi_B, ext_B_val, 1e-8) << "B diverged at step " << step;
    EXPECT_NEAR(bi_C, ext_C_val, 1e-8) << "C diverged at step " << step;

    EXPECT_NEAR(K_EQ * bi_B - bi_C, 0.0, 1e-6) << "Built-in constraint violated at step " << step;
    EXPECT_NEAR(K_EQ * ext_B_val - ext_C_val, 0.0, 1e-6) << "External constraint violated at step " << step;
  }
}

/// @brief Two coupled equilibria using composed models: separate EquilibriumConstraintModels
///        for each ratio constraint plus a MassConservationModel to enforce mass balance.
///
/// System: A → B (driver), B ⇌ C (K1=5), B ⇌ D (K2=3)
/// Constraints: C=K1*B (equilibrium1), D=K2*B (equilibrium2), A+B+C+D=1 (conservation)
/// This makes B, C, D all algebraic; A is the only ODE variable.
/// Steady state: [A]=0, [B]=1/(1+K1+K2), [C]=K1/(1+K1+K2), [D]=K2/(1+K1+K2)
TEST(ExternalModelConstraints, MultiEquilibriumKineticVsComposedConstraints)
{
  constexpr double K1 = 5.0;
  constexpr double K2 = 3.0;
  constexpr double k_drive = 0.1;
  constexpr double k_f1 = 10.0;
  constexpr double k_b1 = k_f1 / K1;
  constexpr double k_f2 = 8.0;
  constexpr double k_b2 = k_f2 / K2;

  double expected_B = 1.0 / (1.0 + K1 + K2);
  double expected_C = K1 * expected_B;
  double expected_D = K2 * expected_B;

  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  auto D = micm::Species("D");
  micm::Phase gas_phase{ "gas", { A, B, C, D } };

  auto system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });

  // ── Kinetic system ──
  micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                             .SetReactants({ A })
                             .SetProducts({ { B, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_drive, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();
  micm::Process rxn_bc = micm::ChemicalReactionBuilder()
                             .SetReactants({ B })
                             .SetProducts({ { C, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_f1, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();
  micm::Process rxn_cb = micm::ChemicalReactionBuilder()
                             .SetReactants({ C })
                             .SetProducts({ { B, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_b1, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();
  micm::Process rxn_bd = micm::ChemicalReactionBuilder()
                             .SetReactants({ B })
                             .SetProducts({ { D, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_f2, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();
  micm::Process rxn_db = micm::ChemicalReactionBuilder()
                             .SetReactants({ D })
                             .SetProducts({ { B, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = k_b2, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();

  auto kin_options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto kin_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(kin_options)
                        .SetSystem(system)
                        .SetReactions({ rxn_ab, rxn_bc, rxn_cb, rxn_bd, rxn_db })
                        .SetReorderState(false)
                        .Build();

  // ── Constraint system: 3 composed external models ──
  EquilibriumConstraintModel eq_bc("B", "C", K1);        // C row: K1*[B]-[C]=0
  EquilibriumConstraintModel eq_bd("B", "D", K2);        // D row: K2*[B]-[D]=0
  MassConservationModel conservation("B", { "A", "B", "C", "D" }, 1.0);  // B row: A+B+C+D-1=0

  auto dae_options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto ext_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(dae_options)
                        .SetSystem(system)
                        .SetReactions({ rxn_ab })
                        .AddExternalModel(eq_bc)
                        .AddExternalModel(eq_bd)
                        .AddExternalModel(conservation)
                        .SetReorderState(false)
                        .Build();

  // ── Solve both ──
  auto state_kin = kin_solver.GetState(1);
  auto state_ext = ext_solver.GetState(1);

  for (auto* s : { &state_kin, &state_ext })
  {
    s->variables_[0][s->variable_map_.at("A")] = 1.0;
    s->variables_[0][s->variable_map_.at("B")] = 0.0;
    s->variables_[0][s->variable_map_.at("C")] = 0.0;
    s->variables_[0][s->variable_map_.at("D")] = 0.0;
    s->conditions_[0].temperature_ = 298.0;
    s->conditions_[0].pressure_ = 101325.0;
  }

  double dt = 1.0;
  for (int step = 0; step < 200; ++step)
  {
    kin_solver.UpdateStateParameters(state_kin);
    ext_solver.UpdateStateParameters(state_ext);

    auto res_kin = kin_solver.Solve(dt, state_kin);
    auto res_ext = ext_solver.Solve(dt, state_ext);

    EXPECT_EQ(res_kin.state_, micm::SolverState::Converged) << "Kinetic failed step " << step;
    EXPECT_EQ(res_ext.state_, micm::SolverState::Converged) << "Constraint failed step " << step;
  }

  double kin_A = state_kin.variables_[0][state_kin.variable_map_.at("A")];
  double kin_B = state_kin.variables_[0][state_kin.variable_map_.at("B")];
  double kin_C = state_kin.variables_[0][state_kin.variable_map_.at("C")];
  double kin_D = state_kin.variables_[0][state_kin.variable_map_.at("D")];

  double ext_A_val = state_ext.variables_[0][state_ext.variable_map_.at("A")];
  double ext_B_val = state_ext.variables_[0][state_ext.variable_map_.at("B")];
  double ext_C_val = state_ext.variables_[0][state_ext.variable_map_.at("C")];
  double ext_D_val = state_ext.variables_[0][state_ext.variable_map_.at("D")];

  // Both near analytical steady state
  EXPECT_NEAR(kin_A, 0.0, 1e-6);
  EXPECT_NEAR(ext_A_val, 0.0, 1e-6);
  EXPECT_NEAR(kin_B, expected_B, 1e-3);
  EXPECT_NEAR(ext_B_val, expected_B, 1e-3);
  EXPECT_NEAR(kin_C, expected_C, 1e-3);
  EXPECT_NEAR(ext_C_val, expected_C, 1e-3);
  EXPECT_NEAR(kin_D, expected_D, 1e-3);
  EXPECT_NEAR(ext_D_val, expected_D, 1e-3);

  // Kinetic and constraint systems agree
  EXPECT_NEAR(kin_B, ext_B_val, 1e-3);
  EXPECT_NEAR(kin_C, ext_C_val, 1e-3);
  EXPECT_NEAR(kin_D, ext_D_val, 1e-3);

  // Mass conservation
  EXPECT_NEAR(kin_A + kin_B + kin_C + kin_D, 1.0, 1e-6);
  EXPECT_NEAR(ext_A_val + ext_B_val + ext_C_val + ext_D_val, 1.0, 1e-3);

  // Equilibrium constraints satisfied
  EXPECT_NEAR(K1 * ext_B_val - ext_C_val, 0.0, 1e-6);
  EXPECT_NEAR(K2 * ext_B_val - ext_D_val, 0.0, 1e-6);
}

/// @brief Verify that process Jacobian elements in algebraic rows not covered by constraints
/// don't cause "Zero element access" exceptions.
///
/// StubAerosolWithSolvent declares process element (A_AQ, S) but the constraint only
/// declares (A_AQ, A_GAS) and (A_AQ, A_AQ). Without the fix, (A_AQ, S) is filtered
/// from the sparsity pattern because A_AQ is algebraic, and the external model's
/// JacobianFunction closure throws when accessing jac[block][A_AQ][S].
TEST(ExternalModelConstraints, ProcessJacobianElementInAlgebraicRowSurvivesFiltering)
{
  auto A_GAS = micm::Species("A_GAS");
  auto S = micm::Species("S");
  micm::Phase gas_phase{ "gas", { A_GAS, S } };

  double total = 1.0;
  double k = 0.1;
  StubAerosolWithSolvent aerosol(k, total);

  auto system = micm::System(micm::SystemParameters{
      .gas_phase_ = gas_phase,
      .external_models_ = { aerosol } });

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

  // This Build() call would throw "Zero element access" without the fix
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(system)
                    .SetReactions({})
                    .AddExternalModel(aerosol)
                    .SetReorderState(false)
                    .Build();

  auto state = solver.GetState(1);

  // Verify A_AQ is algebraic
  auto i_aq = state.variable_map_.at("AEROSOL.A_AQ");
  EXPECT_DOUBLE_EQ(state.upper_left_identity_diagonal_[i_aq], 0.0);

  // Initialize state
  state.variables_[0][state.variable_map_.at("A_GAS")] = 0.9;
  state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")] = 0.1;
  state.variables_[0][state.variable_map_.at("S")] = 1.0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;

  // Solve several steps — verifies no runtime exceptions from Jacobian access
  double dt = 10.0;
  for (int step = 0; step < 10; ++step)
  {
    auto result = solver.Solve(dt, state);
    EXPECT_EQ(result.state_, micm::SolverState::Converged) << "Step " << step;

    double sum = state.variables_[0][state.variable_map_.at("A_GAS")]
               + state.variables_[0][state.variable_map_.at("AEROSOL.A_AQ")];
    EXPECT_NEAR(sum, total, 1e-4) << "Conservation violated at step " << step;
  }
}

// ═══════════════════════════════════════════════════════════════
// Finite-Difference Jacobian Verification for External Models
// ═══════════════════════════════════════════════════════════════

using DenseMatrix = micm::Matrix<double>;
using SparseMatrixFD = micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>;

/// Verify StubAerosolWithConstraints process ForcingFunction/JacobianFunction
TEST(ExternalModelFiniteDifferenceJacobian, ProcessForcingJacobian)
{
  double k = 3.5;
  StubAerosolWithConstraints aerosol(k);

  std::unordered_map<std::string, std::size_t> var_map = { { "A_GAS", 0 }, { "AEROSOL.A_AQ", 1 } };
  std::unordered_map<std::string, std::size_t> param_map;
  const std::size_t num_species = 2;

  auto forcing_fn = aerosol.ForcingFunction<DenseMatrix>(param_map, var_map);
  auto nz_elements = aerosol.NonZeroJacobianElements(var_map);

  auto builder = SparseMatrixFD::Create(num_species).SetNumberOfBlocks(2).InitialValue(0.0);
  for (auto& elem : nz_elements)
    builder = builder.WithElement(elem.first, elem.second);
  SparseMatrixFD analytical_jac{ builder };

  auto jacobian_fn = aerosol.JacobianFunction<DenseMatrix, SparseMatrixFD>(param_map, var_map, analytical_jac);

  DenseMatrix variables(2, num_species, 0.0);
  variables[0][0] = 0.8;
  variables[0][1] = 0.2;
  variables[1][0] = 0.3;
  variables[1][1] = 0.7;

  DenseMatrix params(2, 0, 0.0);

  jacobian_fn(params, variables, analytical_jac);

  auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
  { forcing_fn(params, vars, forcing); };

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

  auto comparison =
      micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

  EXPECT_TRUE(comparison.passed) << "Process Jacobian mismatch: block=" << comparison.worst_block
                                 << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                 << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

  auto sparsity =
      micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

  EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block << " row=" << sparsity.worst_row
                               << " col=" << sparsity.worst_col << " fd_value=" << sparsity.worst_fd;
}

/// Verify StubAerosolWithConstraints constraint ConstraintResidualFunction/ConstraintJacobianFunction
TEST(ExternalModelFiniteDifferenceJacobian, ConstraintResidualJacobian)
{
  double k = 0.1;
  double total = 1.0;
  StubAerosolWithConstraints aerosol(k, total);

  std::unordered_map<std::string, std::size_t> param_map;
  std::unordered_map<std::string, std::size_t> var_map = { { "A_GAS", 0 }, { "AEROSOL.A_AQ", 1 } };
  const std::size_t num_species = 2;

  auto residual_fn = aerosol.ConstraintResidualFunction<DenseMatrix>(param_map, var_map);
  auto nz_elements = aerosol.NonZeroConstraintJacobianElements(var_map);

  auto builder = SparseMatrixFD::Create(num_species).SetNumberOfBlocks(2).InitialValue(0.0);
  for (auto& elem : nz_elements)
    builder = builder.WithElement(elem.first, elem.second);
  SparseMatrixFD analytical_jac{ builder };

  auto jacobian_fn = aerosol.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(param_map, var_map, analytical_jac);

  DenseMatrix variables(2, num_species, 0.0);
  variables[0][0] = 0.6;
  variables[0][1] = 0.4;
  variables[1][0] = 0.2;
  variables[1][1] = 0.8;
  DenseMatrix dummy_params(2, 1, 0.0);

  jacobian_fn(variables, dummy_params, analytical_jac);

  auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
  { residual_fn(vars, dummy_params, forcing); };

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

  auto comparison =
      micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

  EXPECT_TRUE(comparison.passed) << "Constraint Jacobian mismatch: block=" << comparison.worst_block
                                 << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                 << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;
}

/// Verify EquilibriumConstraintModel constraint residual/Jacobian pair
TEST(ExternalModelFiniteDifferenceJacobian, EquilibriumConstraintModelJacobian)
{
  double K_eq = 2.5;
  EquilibriumConstraintModel model("A", "B", K_eq);

  std::unordered_map<std::string, std::size_t> param_map;
  std::unordered_map<std::string, std::size_t> var_map = { { "A", 0 }, { "B", 1 } };
  const std::size_t num_species = 2;

  auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(param_map, var_map);
  auto nz_elements = model.NonZeroConstraintJacobianElements(var_map);

  auto builder = SparseMatrixFD::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);
  for (auto& elem : nz_elements)
    builder = builder.WithElement(elem.first, elem.second);
  SparseMatrixFD analytical_jac{ builder };

  auto jacobian_fn = model.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(param_map, var_map, analytical_jac);

  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;
  variables[0][1] = 5.0;
  DenseMatrix dummy_params(1, 1, 0.0);

  jacobian_fn(variables, dummy_params, analytical_jac);

  auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
  { residual_fn(vars, dummy_params, forcing); };

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

  auto comparison =
      micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

  EXPECT_TRUE(comparison.passed) << "EquilibriumConstraintModel Jacobian mismatch: block=" << comparison.worst_block
                                 << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                 << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

  auto sparsity =
      micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

  EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block << " row=" << sparsity.worst_row
                               << " col=" << sparsity.worst_col << " fd_value=" << sparsity.worst_fd;
}

/// @brief External model constraint with a temperature-dependent K_eq stored as a state parameter
///
/// Enforces: K_eq(T) * [reactant] - [product] = 0
/// where K_eq(T) = K_eq_ref * exp(delta_H / R * (1/T_ref - 1/T))
///
/// This exercises the constraint state parameter pipeline: the model declares a parameter name,
/// provides an update function that computes K_eq from temperature, and the residual/Jacobian
/// functions read K_eq from the state parameter matrix.
class TemperatureDependentEquilibriumModel
{
 public:
  TemperatureDependentEquilibriumModel(
      const std::string& reactant,
      const std::string& product,
      double K_eq_ref,
      double delta_H_over_R,
      double T_ref = 298.15)
      : reactant_(reactant),
        product_(product),
        K_eq_ref_(K_eq_ref),
        delta_H_over_R_(delta_H_over_R),
        T_ref_(T_ref),
        param_name_(product + "_K_eq")
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { product_ };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { reactant_, product_ };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    auto i_r = state_indices.at(reactant_);
    auto i_p = state_indices.at(product_);
    return { { i_p, i_r }, { i_p, i_p } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return { param_name_ };
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>& param_indices) const
  {
    auto i_K = param_indices.at(param_name_);
    double K_ref = K_eq_ref_;
    double dH_R = delta_H_over_R_;
    double T_ref = T_ref_;
    return [=](const std::vector<micm::Conditions>& conditions, DenseMatrixPolicy& params)
    {
      for (std::size_t i = 0; i < conditions.size(); ++i)
      {
        double T = conditions[i].temperature_;
        params[i][i_K] = K_ref * std::exp(dH_R * (1.0 / T_ref - 1.0 / T));
      }
    };
  }

  /// Residual: G = K_eq(T) * [reactant] - [product]
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>& param_indices,
      const std::unordered_map<std::string, std::size_t>& var) const
  {
    auto i_r = var.at(reactant_);
    auto i_p = var.at(product_);
    auto i_K = param_indices.at(param_name_);
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy& params, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
        forcing[i][i_p] = params[i][i_K] * state[i][i_r] - state[i][i_p];
    };
  }

  /// Jacobian: dG/d[reactant] = K_eq(T), dG/d[product] = -1
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>& param_indices,
      const std::unordered_map<std::string, std::size_t>& var,
      const SparseMatrixPolicy&) const
  {
    auto i_r = var.at(reactant_);
    auto i_p = var.at(product_);
    auto i_K = param_indices.at(param_name_);
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& params, SparseMatrixPolicy& jac)
    {
      for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
      {
        jac[i][i_p][i_r] -= params[i][i_K];
        jac[i][i_p][i_p] -= -1.0;
      }
    };
  }

 private:
  std::string reactant_;
  std::string product_;
  double K_eq_ref_;
  double delta_H_over_R_;
  double T_ref_;
  std::string param_name_;
};

/// @brief Verify that external model constraints can use temperature-dependent state parameters
///
/// System: A → B (kinetic), K_eq(T) * [B] - [C] = 0 (algebraic)
/// At T=298.15 K, K_eq = K_eq_ref. At T=350 K, K_eq shifts.
/// The test solves at two temperatures and verifies that [C]/[B] = K_eq(T) at each.
TEST(ExternalModelConstraints, TemperatureDependentConstraintParameter)
{
  auto A = micm::Species("A");
  auto B = micm::Species("B");
  auto C = micm::Species("C");
  micm::Phase gas_phase{ "gas", { A, B, C } };

  const double K_DRIVE = 0.1;
  const double K_EQ_REF = 2.0;
  const double DELTA_H_OVER_R = 3000.0;  // Positive => K_eq increases with T
  const double T_REF = 298.15;

  micm::Process rxn_ab = micm::ChemicalReactionBuilder()
                             .SetReactants({ A })
                             .SetProducts({ { B, 1 } })
                             .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = K_DRIVE, .B_ = 0, .C_ = 0 }))
                             .SetPhase(gas_phase)
                             .Build();

  TemperatureDependentEquilibriumModel eq_model("B", "C", K_EQ_REF, DELTA_H_OVER_R, T_REF);

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ rxn_ab })
                    .AddExternalModel(eq_model)
                    .SetReorderState(false)
                    .Build();

  // Solve at T = 298.15 K  (K_eq = K_EQ_REF = 2.0)
  {
    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.variables_[0][state.variable_map_.at("C")] = 0.0;
    state.conditions_[0].temperature_ = T_REF;
    state.conditions_[0].pressure_ = 101325.0;

    double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged) << "T=298 solve failed at step " << step;
    }

    double B_val = state.variables_[0][state.variable_map_.at("B")];
    double C_val = state.variables_[0][state.variable_map_.at("C")];
    double K_eq_expected = K_EQ_REF;
    EXPECT_GT(B_val, 0.0);
    EXPECT_NEAR(C_val / B_val, K_eq_expected, 1e-4) << "At T=298.15K, [C]/[B] should equal K_eq_ref";
  }

  // Solve at T = 350 K  (K_eq > K_EQ_REF due to positive delta_H)
  {
    auto state = solver.GetState(1);
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.variables_[0][state.variable_map_.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 350.0;
    state.conditions_[0].pressure_ = 101325.0;

    double dt = 1.0;
    for (int step = 0; step < 200; ++step)
    {
      solver.UpdateStateParameters(state);
      auto result = solver.Solve(dt, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged) << "T=350 solve failed at step " << step;
    }

    double B_val = state.variables_[0][state.variable_map_.at("B")];
    double C_val = state.variables_[0][state.variable_map_.at("C")];
    double K_eq_350 = K_EQ_REF * std::exp(DELTA_H_OVER_R * (1.0 / T_REF - 1.0 / 350.0));
    EXPECT_GT(B_val, 0.0);
    EXPECT_GT(K_eq_350, K_EQ_REF) << "K_eq should increase with temperature for positive delta_H";
    EXPECT_NEAR(C_val / B_val, K_eq_350, 1e-4) << "At T=350K, [C]/[B] should equal K_eq(350)";
  }
}
