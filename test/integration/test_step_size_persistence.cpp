// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for RosenbrockSolverParameters::h_persist_: carrying the step-size
// controller's suggestion across public Solve() calls via
// State::solver_step_size_suggestion_.
//
// Without persistence, every Solve() call restarts from h_start (default
// DEFAULT_H_START * time_step) and spends several steps re-growing the step
// size; the controller-cadence experiment measured 9 accepted steps inflating
// to 9000 across 1000 calls. With persistence, a segmented solve should cost
// approximately max(number of segments, single-call steps).

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace micm;

namespace
{
  struct SegmentedResult
  {
    std::uint64_t accepted = 0;
    double final_A = 0.0;
    bool converged = true;
  };

  // First-order decay A -> B with k = 1 s^-1, solved over `total_time` in
  // `segments` equal public Solve() calls.
  SegmentedResult RunDecay(bool h_persist, int segments, double total_time)
  {
    auto A = Species("A");
    auto B = Species("B");
    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B } };

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { B, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0 })
                      .SetPhase(gas_phase)
                      .Build();

    auto options = RosenbrockSolverParameters::FourStageRosenbrockParameters();
    options.h_persist_ = h_persist;
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    SegmentedResult out;
    const double dt = total_time / segments;
    for (int segment = 0; segment < segments; ++segment)
    {
      double advanced = 0.0;
      int guard = 0;
      while (advanced < dt && guard++ < 10000)
      {
        auto result = solver.Solve(dt - advanced, state);
        if (result.state_ != SolverState::Converged)
        {
          out.converged = false;
          return out;
        }
        advanced += result.stats_.final_time_;
        out.accepted += result.stats_.accepted_;
      }
    }
    out.final_A = state.variables_[0][state.variable_map_.at("A")];
    return out;
  }

  // Slaved-equilibrium DAE (B = K_eq * A algebraic) over segmented calls, to
  // exercise persistence together with per-call constraint reinitialization.
  SegmentedResult RunSlavedDae(bool h_persist, int segments, double total_time)
  {
    auto A = Species("A");
    auto B = Species("B");
    auto P = Species("P");
    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ A, B, P } };

    const double k = 1.0;
    const double K_eq = 1.0e-5;

    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ A })
                      .SetProducts({ { P, 1 } })
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = k })
                      .SetPhase(gas_phase)
                      .Build();

    std::vector<Constraint> constraints;
    constraints.push_back(EquilibriumConstraint(
        "eq",
        B,
        std::vector<StoichSpecies>{ { A, 1.0 } },
        std::vector<StoichSpecies>{ { B, 1.0 } },
        VantHoffParam{ .K_HLC_ref_ = K_eq, .delta_H_ = 0.0 }));

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    options.h_persist_ = h_persist;
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
    state.variables_[0][state.variable_map_.at("A")] = 1.0;
    state.variables_[0][state.variable_map_.at("B")] = K_eq;
    state.variables_[0][state.variable_map_.at("P")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    SegmentedResult out;
    const double dt = total_time / segments;
    for (int segment = 0; segment < segments; ++segment)
    {
      double advanced = 0.0;
      int guard = 0;
      while (advanced < dt && guard++ < 10000)
      {
        auto result = solver.Solve(dt - advanced, state);
        if (result.state_ != SolverState::Converged)
        {
          out.converged = false;
          return out;
        }
        advanced += result.stats_.final_time_;
        out.accepted += result.stats_.accepted_;
      }
    }
    out.final_A = state.variables_[0][state.variable_map_.at("A")];
    return out;
  }
}  // namespace

TEST(StepSizePersistence, SegmentedOdeCostsNearSingleCall)
{
  constexpr int kSegments = 20;
  constexpr double kTotalTime = 10.0;

  auto single = RunDecay(false, 1, kTotalTime);
  auto segmented_off = RunDecay(false, kSegments, kTotalTime);
  auto segmented_on = RunDecay(true, kSegments, kTotalTime);

  ASSERT_TRUE(single.converged);
  ASSERT_TRUE(segmented_off.converged);
  ASSERT_TRUE(segmented_on.converged);

  // Persistence must recover most of the segmentation overhead: bounded by the
  // single-call cost plus roughly one step per segment (plus slack for the
  // first-segment ramp), and strictly cheaper than the restarting behavior.
  EXPECT_LT(segmented_on.accepted, segmented_off.accepted);
  EXPECT_LE(segmented_on.accepted, single.accepted + kSegments + 5);

  // All three must agree with the analytic solution exp(-t).
  const double exact = std::exp(-kTotalTime);
  EXPECT_NEAR(single.final_A, exact, 1.0e-5 * exact + 1.0e-12);
  EXPECT_NEAR(segmented_off.final_A, exact, 1.0e-5 * exact + 1.0e-12);
  EXPECT_NEAR(segmented_on.final_A, exact, 1.0e-5 * exact + 1.0e-12);
}

TEST(StepSizePersistence, FlagOffMatchesRestartingBaseline)
{
  // With the flag off (default), a fresh state must behave identically on
  // repeated runs: persistence state must not leak into default behavior.
  auto first = RunDecay(false, 5, 10.0);
  auto second = RunDecay(false, 5, 10.0);
  ASSERT_TRUE(first.converged);
  EXPECT_EQ(first.accepted, second.accepted);
  EXPECT_EQ(first.final_A, second.final_A);
}

TEST(StepSizePersistence, SegmentedDaeConvergesAndAgrees)
{
  constexpr int kSegments = 20;
  constexpr double kTotalTime = 10.0;

  auto segmented_off = RunSlavedDae(false, kSegments, kTotalTime);
  auto segmented_on = RunSlavedDae(true, kSegments, kTotalTime);

  ASSERT_TRUE(segmented_off.converged);
  ASSERT_TRUE(segmented_on.converged);

  // Same trajectory to integration accuracy, at reduced (or equal) cost;
  // per-call constraint reinitialization must still succeed when the carried
  // step size skips the ramp.
  EXPECT_LE(segmented_on.accepted, segmented_off.accepted);
  const double rel_diff = std::abs(segmented_on.final_A - segmented_off.final_A) /
                          (std::abs(segmented_off.final_A) + 1.0e-30);
  EXPECT_LT(rel_diff, 1.0e-4);
}
