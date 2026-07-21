// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for the step-acceptance error-norm policy (State::cellwise_error_norm_)
// and the constraint-initialization pivot diagnostic
// (SolverStats::constraint_init_min_pivot_ratio_).
//
// Norm policy: with a global WRMS, one active cell's error is diluted by the
// quiescent cells sharing the batch (sqrt(N) fewer effective digits at 100
// stationary cells); the cellwise-max norm must hold the active cell's
// accuracy batch-invariant.
//
// Pivot diagnostic: a converged Newton correction cannot certify forward error
// once the algebraic block is numerically singular; the solver must surface a
// small algebraic-pivot ratio so callers can flag such constraints.

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace micm;

namespace
{
  struct DecayResult
  {
    std::uint64_t accepted = 0;
    double active_cell_error = 0.0;
    bool converged = true;
  };

  // First-order decay A -> B in `cells` cells; only the last cell is active
  // (A = 1), all others are exactly stationary (A = 0). The active cell's
  // relative error against exp(-t) is the metric.
  DecayResult RunDilution(bool cellwise_norm, std::size_t cells, double total_time)
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
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(cells);
    state.cellwise_error_norm_ = cellwise_norm;
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-14));
    const auto a = state.variable_map_.at("A");
    const auto b = state.variable_map_.at("B");
    for (std::size_t cell = 0; cell < cells; ++cell)
    {
      state.variables_[cell][a] = 0.0;
      state.variables_[cell][b] = 0.0;
      state.conditions_[cell].temperature_ = 298.0;
      state.conditions_[cell].pressure_ = 101325.0;
    }
    state.variables_[cells - 1][a] = 1.0;
    solver.UpdateStateParameters(state);

    DecayResult out;
    double advanced = 0.0;
    int guard = 0;
    while (advanced < total_time && guard++ < 10000)
    {
      auto result = solver.Solve(total_time - advanced, state);
      if (result.state_ != SolverState::Converged)
      {
        out.converged = false;
        return out;
      }
      advanced += result.stats_.final_time_;
      out.accepted += result.stats_.accepted_;
    }
    const double exact = std::exp(-total_time);
    out.active_cell_error = std::abs(state.variables_[cells - 1][a] - exact) / exact;
    return out;
  }

  /// Constraint-only external model with a nearly singular algebraic row:
  /// G = epsilon * Z - X. dG/dZ = epsilon, so the algebraic U pivot is
  /// epsilon while the differential (identity) pivots are 1.
  class NearSingularConstraintModel
  {
   public:
    explicit NearSingularConstraintModel(double epsilon)
        : epsilon_(epsilon)
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto z = s.at("Z");
      return { { z, s.at("X") }, { z, z } };
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x = s.at("X");
      const auto z = s.at("Z");
      const double epsilon = epsilon_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][z] = epsilon * state[cell][z] - state[cell][x];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto x = s.at("X");
      const auto z = s.at("Z");
      const double epsilon = epsilon_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][z][x] -= -1.0;
          jacobian[block][z][z] -= epsilon;
        }
      };
    }

   private:
    double epsilon_;
  };

  double PivotRatioFor(double epsilon)
  {
    auto X = Species("X");
    auto Z = Species("Z");
    Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ X, Z } };
    Process rxn = ChemicalReactionBuilder()
                      .SetReactants({ X })
                      .SetProducts({})
                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0 })
                      .SetPhase(gas_phase)
                      .Build();
    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions({ rxn })
                      .AddExternalModel(NearSingularConstraintModel(epsilon))
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
    // Off the manifold so initialization must factor at least once.
    state.variables_[0][state.variable_map_.at("X")] = 1.0;
    state.variables_[0][state.variable_map_.at("Z")] = 0.0;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(1.0e-3, state);
    EXPECT_EQ(result.state_, SolverState::Converged);
    return result.stats_.constraint_init_min_pivot_ratio_;
  }
}  // namespace

TEST(CellwiseErrorNorm, ActiveCellAccuracyIsBatchInvariant)
{
  constexpr double kTotalTime = 1.0;
  const auto single = RunDilution(false, 1, kTotalTime);
  const auto global_100 = RunDilution(false, 100, kTotalTime);
  const auto cellwise_100 = RunDilution(true, 100, kTotalTime);

  ASSERT_TRUE(single.converged);
  ASSERT_TRUE(global_100.converged);
  ASSERT_TRUE(cellwise_100.converged);

  // The global norm dilutes the active cell's error control by ~sqrt(100);
  // the cellwise norm must keep the active cell's accuracy near the
  // single-cell result (within a small factor) and take at least as many
  // steps as the diluted global run.
  EXPECT_GT(global_100.active_cell_error, 3.0 * single.active_cell_error);
  EXPECT_LT(cellwise_100.active_cell_error, 3.0 * single.active_cell_error);
  EXPECT_GE(cellwise_100.accepted, global_100.accepted);
}

TEST(ConstraintInitPivotDiagnostic, FlagsNearSingularAlgebraicBlock)
{
  const double healthy = PivotRatioFor(1.0);
  const double singular = PivotRatioFor(1.0e-12);

  // Healthy constraint: algebraic pivot comparable to the identity rows.
  EXPECT_GT(healthy, 1.0e-2);
  // Nearly singular constraint: ratio collapses with epsilon.
  EXPECT_LT(singular, 1.0e-10);
}
