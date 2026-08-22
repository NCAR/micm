// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

// CPU/Kokkos numerical agreement.
//
// Every other test in this directory checks the Kokkos backend against an analytical solution
// or an invariant. None of them check it against the CPU backend, so a Kokkos-only defect that
// stays inside the solver's error tolerance -- a mis-indexed Jacobian entry, a stale device
// buffer, a padding lane leaking into a reduction -- can hide behind a loose tolerance.
//
// This test runs the same mechanism, the same initial conditions and the same solver
// parameters through both backends and compares the trajectories step by step.
//
// The CPU side is deliberately built with CpuSolverBuilderInPlace over
// LuDecompositionMozartInPlace, because micm::KokkosSolverBuilder expands to
//   SolverBuilder<Params, KokkosDenseMatrix<Real,L>, KokkosSparseMatrix<Real, vector ordering>,
//                 ProcessSet<...>, LuDecompositionMozartInPlace<...>,
//                 LinearSolverInPlace<...>, State<...>>
// (see include/micm/kokkos/solver/kokkos_solver_builder.hpp) and CpuSolverBuilderInPlace
// expands to the identical shape over VectorMatrix / SparseMatrix. Matching every policy is
// the point: the two solvers then perform the same operations in the same order, so any
// difference is attributable to the Kokkos matrix implementation rather than to a different
// algorithm.

#include "../../precision_matchers.hpp"

#include <micm/Kokkos.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{
  constexpr micm::Index NUM_CELLS = 5;

  template<micm::Index L>
  using CpuSparse = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>;

  template<micm::Index L>
  using CpuBuilder = micm::CpuSolverBuilderInPlace<
      micm::RosenbrockSolverParameters,
      micm::VectorMatrix<micm::Real, L>,
      CpuSparse<L>,
      micm::LuDecompositionMozartInPlace<CpuSparse<L>>>;

  template<micm::Index L>
  using KokkosBuilder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters, L>;

  /// @brief One step's worth of results, keyed by name so the two backends may order the
  ///        state differently without the comparison caring.
  struct Snapshot
  {
    micm::SolverState solver_state_{ micm::SolverState::NotYetCalled };
    micm::Real final_time_{};
    std::map<std::string, std::vector<micm::Real>> variables_;
    std::vector<std::vector<micm::Real>> rate_constants_;
  };

  template<class StateType>
  Snapshot Capture(StateType& state, const micm::SolverResult& result)
  {
    Snapshot snapshot;
    snapshot.solver_state_ = result.state_;
    snapshot.final_time_ = result.stats_.final_time_;
    for (const auto& [name, index] : state.variable_map_)
    {
      std::vector<micm::Real> per_cell;
      per_cell.reserve(state.variables_.NumRows());
      for (micm::Index cell = 0; cell < state.variables_.NumRows(); ++cell)
      {
        per_cell.push_back(state.variables_[cell][index]);
      }
      snapshot.variables_.emplace(name, std::move(per_cell));
    }
    snapshot.rate_constants_.reserve(state.rate_constants_.NumRows());
    for (micm::Index cell = 0; cell < state.rate_constants_.NumRows(); ++cell)
    {
      std::vector<micm::Real> per_reaction;
      per_reaction.reserve(state.rate_constants_.NumColumns());
      for (micm::Index i = 0; i < state.rate_constants_.NumColumns(); ++i)
      {
        per_reaction.push_back(state.rate_constants_[cell][i]);
      }
      snapshot.rate_constants_.push_back(std::move(per_reaction));
    }
    return snapshot;
  }

  /// @brief Largest relative difference seen so far, reported at the end of each comparison so
  ///        the tolerance below can be justified from measurement rather than guessed.
  double g_max_relative_difference = 0.0;

  void ExpectAgree(micm::Real cpu, micm::Real kokkos, const std::string& what)
  {
    // The two backends should be performing identical arithmetic, so the bar is far tighter
    // than the solver's own relative tolerance (1e-8 in double). It is not set to zero
    // because a device kernel is entitled to reassociate a reduction.
    //
    // Measured with the Serial backend in a double build: the largest relative difference over
    // all three mechanisms and L = 1, 2, 4, 8 is 4.8e-15, and L=1 plus the constrained case are
    // bit-identical. 1e-10 therefore leaves four orders of magnitude of headroom for a CUDA
    // backend to reassociate more aggressively, while still being 100x tighter than the
    // solver tolerance that any genuine backend defect would have to hide inside.
    constexpr double rel_tol = std::is_same_v<micm::Real, double> ? 1.0e-10 : 1.0e-4;
    constexpr double abs_floor = std::is_same_v<micm::Real, double> ? 1.0e-18 : 1.0e-12;

    const auto a = static_cast<double>(cpu);
    const auto b = static_cast<double>(kokkos);
    if (std::abs(a) > abs_floor)
    {
      g_max_relative_difference = std::max(g_max_relative_difference, std::abs(a - b) / std::abs(a));
    }
    EXPECT_NEAR(b, a, abs_floor + rel_tol * std::abs(a)) << what;
  }

  void ExpectTrajectoriesAgree(
      const std::vector<Snapshot>& cpu,
      const std::vector<Snapshot>& kokkos,
      const std::string& label)
  {
    ASSERT_EQ(cpu.size(), kokkos.size()) << label << ": different number of steps captured";
    for (std::size_t step = 0; step < cpu.size(); ++step)
    {
      const auto& c = cpu[step];
      const auto& k = kokkos[step];
      std::string where = label;
      where.append(" step ").append(std::to_string(step));

      EXPECT_EQ(c.solver_state_, micm::SolverState::Converged)
          << where << ": CPU solver reported " << micm::SolverStateToString(c.solver_state_);
      EXPECT_EQ(k.solver_state_, micm::SolverState::Converged)
          << where << ": Kokkos solver reported " << micm::SolverStateToString(k.solver_state_);

      // Both were asked to advance the same interval, so the integrated time must match even
      // if the internal step-size histories differ.
      ExpectAgree(c.final_time_, k.final_time_, where + ": final_time_");

      ASSERT_EQ(c.variables_.size(), k.variables_.size()) << where << ": different species counts";
      for (const auto& [name, cpu_cells] : c.variables_)
      {
        const auto it = k.variables_.find(name);
        ASSERT_NE(it, k.variables_.end()) << where << ": Kokkos state has no species " << name;
        ASSERT_EQ(cpu_cells.size(), it->second.size()) << where << ": cell count differs for " << name;
        for (std::size_t cell = 0; cell < cpu_cells.size(); ++cell)
        {
          std::string cell_label = where;
          cell_label.append(": ").append(name).append(" cell ").append(std::to_string(cell));
          ExpectAgree(cpu_cells[cell], it->second[cell], cell_label);
        }
      }

      ASSERT_EQ(c.rate_constants_.size(), k.rate_constants_.size()) << where << ": rate-constant cell count differs";
      for (std::size_t cell = 0; cell < c.rate_constants_.size(); ++cell)
      {
        ASSERT_EQ(c.rate_constants_[cell].size(), k.rate_constants_[cell].size())
            << where << ": rate-constant count differs in cell " << cell;
        for (std::size_t i = 0; i < c.rate_constants_[cell].size(); ++i)
        {
          std::string rate_label = where;
          rate_label.append(": rate constant ").append(std::to_string(i)).append(" cell ").append(std::to_string(cell));
          ExpectAgree(c.rate_constants_[cell][i], k.rate_constants_[cell][i], rate_label);
        }
      }
    }
  }

  /// @brief Assert the trajectory actually moved, so "both backends did nothing" cannot pass.
  void ExpectEvolved(const std::vector<Snapshot>& trajectory, const std::string& species, const std::string& label)
  {
    ASSERT_FALSE(trajectory.empty()) << label << ": no steps captured";
    const auto& first = trajectory.front().variables_.at(species);
    const auto& last = trajectory.back().variables_.at(species);
    double largest_change = 0.0;
    for (std::size_t cell = 0; cell < first.size(); ++cell)
    {
      largest_change = std::max(largest_change, std::abs(static_cast<double>(last[cell] - first[cell])));
    }
    EXPECT_GT(largest_change, 0.0) << label << ": " << species << " never changed -- the comparison is vacuous";
  }

  // ---------------------------------------------------------------------------
  // Mechanism 1: Chapman
  // ---------------------------------------------------------------------------

  template<class BuilderPolicy>
  std::vector<Snapshot> RunChapman(BuilderPolicy builder)
  {
    auto o = micm::Species("O");
    auto o1d = micm::Species("O1D");
    auto o2 = micm::Species("O2");
    auto o3 = micm::Species("O3");
    auto m = micm::Species("M");
    auto n2 = micm::Species("N2");

    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ o, o1d, o2, o3, m, n2 } };

    micm::Process r1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o1d, n2 })
                           .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(n2, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 })
                           .SetPhase(gas_phase)
                           .Build();
    micm::Process r2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o, o3 })
                           .SetProducts({ micm::StoichSpecies(o2, 2) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 })
                           .SetPhase(gas_phase)
                           .Build();
    micm::Process r3 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o, o2, m })
                           .SetProducts({ micm::StoichSpecies(o3, 1), micm::StoichSpecies(m, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 })
                           .SetPhase(gas_phase)
                           .Build();
    micm::Process photo_1 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o2 })
                                .SetProducts({ micm::StoichSpecies(o, 2) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO2" })
                                .SetPhase(gas_phase)
                                .Build();
    micm::Process photo_2 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o3 })
                                .SetProducts({ micm::StoichSpecies(o1d, 1), micm::StoichSpecies(o2, 1) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3a" })
                                .SetPhase(gas_phase)
                                .Build();
    micm::Process photo_3 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o3 })
                                .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(o2, 1) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3b" })
                                .SetPhase(gas_phase)
                                .Build();

    auto solver = builder.SetSystem(micm::System(gas_phase))
                      .SetReactions({ r1, r2, r3, photo_1, photo_2, photo_3 })
                      .SetIgnoreUnusedSpecies(true)
                      .Build();
    auto state = solver.GetState(NUM_CELLS);

    // Distinct conditions per cell, so cross-cell vectorisation and the padding tail are both
    // exercised rather than every lane carrying the same numbers.
    for (micm::Index cell = 0; cell < NUM_CELLS; ++cell)
    {
      const micm::Real scale = 1.0 + 0.25 * static_cast<micm::Real>(cell);
      state.variables_[cell][state.variable_map_.at("O")] = 1.0e-9 * scale;
      state.variables_[cell][state.variable_map_.at("O1D")] = 1.0e-14 * scale;
      state.variables_[cell][state.variable_map_.at("O2")] = 5.0e18 * scale;
      state.variables_[cell][state.variable_map_.at("O3")] = 1.0e12 * scale;
      state.variables_[cell][state.variable_map_.at("M")] = 2.5e19 * scale;
      state.variables_[cell][state.variable_map_.at("N2")] = 2.0e19 * scale;
      state.conditions_[cell].temperature_ = 235.0 + 12.0 * static_cast<micm::Real>(cell);
      state.conditions_[cell].pressure_ = 90000.0 + 2000.0 * static_cast<micm::Real>(cell);
      state.conditions_[cell].air_density_ = 2.5e19 * scale;
    }
    std::vector<micm::Real> j_o2(NUM_CELLS);
    std::vector<micm::Real> j_o3a(NUM_CELLS);
    std::vector<micm::Real> j_o3b(NUM_CELLS);
    for (micm::Index cell = 0; cell < NUM_CELLS; ++cell)
    {
      j_o2[cell] = 1.0e-10 * (1.0 + 0.1 * static_cast<micm::Real>(cell));
      j_o3a[cell] = 1.0e-4 * (1.0 + 0.1 * static_cast<micm::Real>(cell));
      j_o3b[cell] = 5.0e-4 * (1.0 + 0.1 * static_cast<micm::Real>(cell));
    }
    state.SetCustomRateParameter("jO2", j_o2);
    state.SetCustomRateParameter("jO3a", j_o3a);
    state.SetCustomRateParameter("jO3b", j_o3b);

    state.variables_.CopyToDevice();
    state.conditions_.CopyToDevice();
    state.custom_rate_parameters_.CopyToDevice();
    solver.UpdateStateParameters(state);

    std::vector<Snapshot> trajectory;
    for (micm::Index step = 0; step < 10; ++step)
    {
      auto result = solver.Solve(30.0, state);
      state.variables_.CopyToHost();
      state.rate_constants_.CopyToHost();
      trajectory.push_back(Capture(state, result));
    }
    return trajectory;
  }

  // ---------------------------------------------------------------------------
  // Mechanism 2: Robertson (stiff)
  // ---------------------------------------------------------------------------

  template<class BuilderPolicy>
  std::vector<Snapshot> RunRobertson(BuilderPolicy builder)
  {
    auto a = micm::Species("A");
    auto b = micm::Species("B");
    auto c = micm::Species("C");
    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

    micm::Process r1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ a })
                           .SetProducts({ micm::StoichSpecies(b, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r1" })
                           .SetPhase(gas_phase)
                           .Build();
    micm::Process r2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, b })
                           .SetProducts({ micm::StoichSpecies(b, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" })
                           .SetPhase(gas_phase)
                           .Build();
    micm::Process r3 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, c })
                           .SetProducts({ micm::StoichSpecies(a, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r3" })
                           .SetPhase(gas_phase)
                           .Build();

    auto solver =
        builder.SetSystem(micm::System(gas_phase)).SetReactions({ r1, r2, r3 }).SetIgnoreUnusedSpecies(true).Build();
    auto state = solver.GetState(NUM_CELLS);
    state.SetRelativeTolerance(micm::Real{ std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-5 });

    for (micm::Index cell = 0; cell < NUM_CELLS; ++cell)
    {
      state.variables_[cell][state.variable_map_.at("A")] = 1.0 - 0.1 * static_cast<micm::Real>(cell);
      state.variables_[cell][state.variable_map_.at("B")] = 0.0;
      state.variables_[cell][state.variable_map_.at("C")] = 0.1 * static_cast<micm::Real>(cell);
      state.conditions_[cell].temperature_ = 280.0 + 5.0 * static_cast<micm::Real>(cell);
      state.conditions_[cell].pressure_ = 101325.0;
      state.conditions_[cell].air_density_ = 2.5e19;
    }
    state.SetCustomRateParameter("r1", std::vector<micm::Real>(NUM_CELLS, 0.04));
    state.SetCustomRateParameter("r2", std::vector<micm::Real>(NUM_CELLS, 3.0e7));
    state.SetCustomRateParameter("r3", std::vector<micm::Real>(NUM_CELLS, 1.0e4));

    state.variables_.CopyToDevice();
    state.conditions_.CopyToDevice();
    state.custom_rate_parameters_.CopyToDevice();
    solver.UpdateStateParameters(state);

    std::vector<Snapshot> trajectory;
    micm::Real time_step = 1.0e-3;
    for (micm::Index step = 0; step < 10; ++step)
    {
      auto result = solver.Solve(time_step, state);
      state.variables_.CopyToHost();
      state.rate_constants_.CopyToHost();
      trajectory.push_back(Capture(state, result));
      time_step *= 4.0;
    }
    return trajectory;
  }

  // ---------------------------------------------------------------------------
  // Mechanism 3: DAE with a linear constraint, so ConstraintSet's device path is compared
  // ---------------------------------------------------------------------------

  template<class BuilderPolicy>
  std::vector<Snapshot> RunConstrained(BuilderPolicy builder)
  {
    using DenseMatrix = typename BuilderPolicy::DenseMatrixPolicyType;
    using SparseMatrix = typename BuilderPolicy::SparseMatrixPolicyType;

    auto a = micm::Species("A");
    auto b = micm::Species("B");
    auto c = micm::Species("C");
    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

    micm::Process rxn = micm::ChemicalReactionBuilder()
                            .SetReactants({ a })
                            .SetProducts({ micm::StoichSpecies(b, 1) })
                            .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "k_ab" })
                            .SetPhase(gas_phase)
                            .Build();

    const micm::Real total = 1.0e-6;
    std::vector<micm::Constraint<DenseMatrix, SparseMatrix>> constraints;
    constraints.emplace_back(micm::LinearConstraint<DenseMatrix, SparseMatrix>(
        "mass_conservation", c, { { a, 1.0 }, { b, 1.0 }, { c, 1.0 } }, total));

    auto solver = builder.SetSystem(micm::System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetIgnoreUnusedSpecies(true)
                      .Build();
    auto state = solver.GetState(NUM_CELLS);

    for (micm::Index cell = 0; cell < NUM_CELLS; ++cell)
    {
      const micm::Real a0 = 0.9e-6 - 0.1e-6 * static_cast<micm::Real>(cell);
      state.variables_[cell][state.variable_map_.at("A")] = a0;
      state.variables_[cell][state.variable_map_.at("B")] = 0.0;
      state.variables_[cell][state.variable_map_.at("C")] = total - a0;
      state.conditions_[cell].temperature_ = 290.0 + 3.0 * static_cast<micm::Real>(cell);
      state.conditions_[cell].pressure_ = 101325.0;
      state.conditions_[cell].air_density_ = 2.5e19;
    }
    state.SetCustomRateParameter("k_ab", std::vector<micm::Real>(NUM_CELLS, 0.5));

    state.variables_.CopyToDevice();
    state.conditions_.CopyToDevice();
    state.custom_rate_parameters_.CopyToDevice();
    solver.UpdateStateParameters(state);

    std::vector<Snapshot> trajectory;
    for (micm::Index step = 0; step < 10; ++step)
    {
      auto result = solver.Solve(0.5, state);
      state.variables_.CopyToHost();
      state.rate_constants_.CopyToHost();
      trajectory.push_back(Capture(state, result));
    }
    return trajectory;
  }

  template<micm::Index L, class Runner>
  void CompareBackends(Runner run, const std::string& label, const std::string& evolving_species)
  {
    g_max_relative_difference = 0.0;
    auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
    auto cpu = run(CpuBuilder<L>(options));
    auto kokkos = run(KokkosBuilder<L>(options));
    std::string tagged = label;
    tagged.append(" L=").append(std::to_string(L));
    ExpectEvolved(cpu, evolving_species, tagged + " (CPU)");
    ExpectEvolved(kokkos, evolving_species, tagged + " (Kokkos)");
    ExpectTrajectoriesAgree(cpu, kokkos, tagged);
    std::cout << "[ AGREEMENT ] " << tagged << " max relative difference = " << g_max_relative_difference << std::endl;
  }
}  // namespace

TEST(KokkosCpuAgreement, Chapman)
{
  CompareBackends<1>([](auto builder) { return RunChapman(std::move(builder)); }, "Chapman", "O3");
  CompareBackends<2>([](auto builder) { return RunChapman(std::move(builder)); }, "Chapman", "O3");
  CompareBackends<4>([](auto builder) { return RunChapman(std::move(builder)); }, "Chapman", "O3");
  CompareBackends<8>([](auto builder) { return RunChapman(std::move(builder)); }, "Chapman", "O3");
}

TEST(KokkosCpuAgreement, RobertsonStiff)
{
  CompareBackends<1>([](auto builder) { return RunRobertson(std::move(builder)); }, "Robertson", "C");
  CompareBackends<2>([](auto builder) { return RunRobertson(std::move(builder)); }, "Robertson", "C");
  CompareBackends<4>([](auto builder) { return RunRobertson(std::move(builder)); }, "Robertson", "C");
  CompareBackends<8>([](auto builder) { return RunRobertson(std::move(builder)); }, "Robertson", "C");
}

TEST(KokkosCpuAgreement, LinearConstraintDae)
{
  CompareBackends<1>([](auto builder) { return RunConstrained(std::move(builder)); }, "LinearConstraint", "B");
  CompareBackends<2>([](auto builder) { return RunConstrained(std::move(builder)); }, "LinearConstraint", "B");
  CompareBackends<4>([](auto builder) { return RunConstrained(std::move(builder)); }, "LinearConstraint", "B");
  CompareBackends<8>([](auto builder) { return RunConstrained(std::move(builder)); }, "LinearConstraint", "B");
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
