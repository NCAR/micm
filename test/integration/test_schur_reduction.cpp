// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for RosenbrockSolverParameters::schur_reduction_: DAE stage solves on
// the Schur complement of the differential block. The reduction is exact
// linear algebra, so with identical (fixed) step sequences the reduced and
// unreduced solutions must agree to roundoff; with adaptive steps the
// trajectories must agree at the integration tolerance.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using namespace micm;

namespace
{
  // N equilibrium triples (A <-> B fast, B -> C slow) with equilibrium +
  // conservation constraints: two coupled algebraic variables per triple.
  struct TripleResult
  {
    std::vector<double> values;
    std::uint64_t accepted = 0;
    bool converged = true;
  };

  TripleResult RunTriples(int n_triples, std::size_t cells, bool schur, bool fixed_step)
  {
    constexpr double Keq = 1.0, Ctot = 1.0, ks = 1.0, S = 1.0e6;
    std::vector<PhaseSpecies> species;
    for (int i = 0; i < n_triples; ++i)
    {
      species.push_back(Species("A" + std::to_string(i)));
      species.push_back(Species("B" + std::to_string(i)));
      species.push_back(Species("C" + std::to_string(i)));
    }
    Phase gas{ "gas", species };
    std::vector<Process> processes;
    std::vector<Constraint> constraints;
    for (int i = 0; i < n_triples; ++i)
    {
      auto A = Species("A" + std::to_string(i));
      auto B = Species("B" + std::to_string(i));
      auto C = Species("C" + std::to_string(i));
      processes.push_back(ChemicalReactionBuilder()
                              .SetReactants({ B })
                              .SetProducts({ { C, 1 } })
                              .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = ks })
                              .SetPhase(gas)
                              .Build());
      constraints.push_back(EquilibriumConstraint(
          "eq" + std::to_string(i),
          B,
          std::vector<StoichSpecies>{ { A, 1.0 } },
          std::vector<StoichSpecies>{ { B, 1.0 } },
          VantHoffParam{ .K_HLC_ref_ = Keq, .delta_H_ = 0.0 }));
      constraints.push_back(
          LinearConstraint("mass" + std::to_string(i), A, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, Ctot));
    }
    (void)S;

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    options.schur_reduction_ = schur;
    if (fixed_step)
    {
      options.h_min_ = 0.125;
      options.h_max_ = 0.125;
      options.h_start_ = 0.125;
      options.max_number_of_steps_ = 100000;
    }
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas))
                      .SetReactions(processes)
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(cells);
    state.SetRelativeTolerance(fixed_step ? 1.0e6 : 1.0e-8);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, fixed_step ? 1.0e6 : 1.0e-14));
    for (std::size_t cell = 0; cell < cells; ++cell)
    {
      for (int i = 0; i < n_triples; ++i)
      {
        // Slightly cell-dependent manifold start so multi-cell states differ.
        const double total = Ctot * (1.0 + 0.1 * cell);
        state.variables_[cell][state.variable_map_.at("A" + std::to_string(i))] = total / (1.0 + Keq);
        state.variables_[cell][state.variable_map_.at("B" + std::to_string(i))] = total * Keq / (1.0 + Keq);
        state.variables_[cell][state.variable_map_.at("C" + std::to_string(i))] = 0.0;
      }
      state.conditions_[cell].temperature_ = 298.0;
      state.conditions_[cell].pressure_ = 101325.0;
    }
    solver.UpdateStateParameters(state);

    TripleResult out;
    double advanced = 0.0;
    int guard = 0;
    while (advanced < 2.0 && guard++ < 10000)
    {
      auto result = solver.Solve(2.0 - advanced, state);
      if (result.state_ != SolverState::Converged && result.state_ != SolverState::ConvergenceExceededMaxSteps)
      {
        out.converged = false;
        return out;
      }
      out.accepted += result.stats_.accepted_;
      advanced += result.stats_.final_time_;
    }
    for (std::size_t cell = 0; cell < cells; ++cell)
      for (std::size_t v = 0; v < state.state_size_; ++v)
        out.values.push_back(state.variables_[cell][v]);
    return out;
  }

  // Robertson QSSA (single coupled algebraic row) via the external-model path.
  double RunRobertsonQssa(bool schur)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto P = Species("P");
    Phase gas{ "gas", std::vector<PhaseSpecies>{ A, B, P } };
    const double K_eq = 1.0e-5;
    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { P, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0 })
                      .SetPhase(gas)
                      .Build();
    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq",
        B,
        std::vector<StoichSpecies>{ { A, 1.0 } },
        std::vector<StoichSpecies>{ { B, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = 0.0 }));
    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    options.schur_reduction_ = schur;
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas))
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
    while (advanced < 5.0 && guard++ < 10000)
    {
      auto result = solver.Solve(5.0 - advanced, state);
      EXPECT_EQ(result.state_, SolverState::Converged);
      if (result.state_ != SolverState::Converged)
        break;
      advanced += result.stats_.final_time_;
    }
    return state.variables_[0][state.variable_map_.at("A")];
  }
}  // namespace

TEST(SchurReduction, FixedStepSolutionsMatchToRoundoff)
{
  // Identical forced step sequences: the reduction is exact linear algebra,
  // so per-step solutions differ only by factorization roundoff.
  auto full = RunTriples(3, 1, false, true);
  auto reduced = RunTriples(3, 1, true, true);
  ASSERT_TRUE(full.converged);
  ASSERT_TRUE(reduced.converged);
  ASSERT_EQ(full.accepted, reduced.accepted);
  ASSERT_EQ(full.values.size(), reduced.values.size());
  for (std::size_t i = 0; i < full.values.size(); ++i)
  {
    EXPECT_NEAR(full.values[i], reduced.values[i], 1.0e-12 * std::max(1.0, std::abs(full.values[i])));
  }
}

TEST(SchurReduction, AdaptiveTrajectoriesAgreeAtTolerance)
{
  auto full = RunTriples(4, 1, false, false);
  auto reduced = RunTriples(4, 1, true, false);
  ASSERT_TRUE(full.converged);
  ASSERT_TRUE(reduced.converged);
  for (std::size_t i = 0; i < full.values.size(); ++i)
  {
    EXPECT_NEAR(full.values[i], reduced.values[i], 1.0e-7 * std::max(1.0e-3, std::abs(full.values[i])));
  }
}

TEST(SchurReduction, MultiCellStatesAgree)
{
  auto full = RunTriples(2, 3, false, false);
  auto reduced = RunTriples(2, 3, true, false);
  ASSERT_TRUE(full.converged);
  ASSERT_TRUE(reduced.converged);
  for (std::size_t i = 0; i < full.values.size(); ++i)
  {
    EXPECT_NEAR(full.values[i], reduced.values[i], 1.0e-7 * std::max(1.0e-3, std::abs(full.values[i])));
  }
}

TEST(SchurReduction, ExternalModelConstraintAgrees)
{
  const double full = RunRobertsonQssa(false);
  const double reduced = RunRobertsonQssa(true);
  const double exact = std::exp(-5.0);
  EXPECT_NEAR(reduced, full, 1.0e-8 * full);
  EXPECT_NEAR(reduced, exact, 1.0e-6 * exact);
}
