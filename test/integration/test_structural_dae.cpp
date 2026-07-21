// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Structural and adversarial DAE coverage:
//
//   INDEX-2 GUARD: a constraint row that does not involve its own algebraic
//   variable (dG/dz = 0) makes the algebraic block singular — the formulation
//   is not index 1. The solver must fail with a named SolverState and roll the
//   caller's state back, never return a silently wrong answer.
//
//   RANDOMIZED MECHANISM SWEEP: seeded random mass-action mechanisms with
//   randomly chosen QSSA species (each guaranteed a strong linear loss so the
//   row is index 1), integrated as DAEs. Contract: the solver converges, the
//   final state satisfies the constraint residuals relatively, and every value
//   is finite and nonnegative. Failures print the seed for reproduction.

#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <functional>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace micm;

namespace
{
  // ---------------------------------------------------------------------------
  // Index-2 guard: G = X - 0.5 with Z algebraic (dG/dZ = 0).
  // ---------------------------------------------------------------------------
  class NotIndexOneModel
  {
   public:
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
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][z] = state[cell][x] - 0.5;  // no dependence on Z
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
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][z][x] -= 1.0;
          jacobian[block][z][z] -= 0.0;  // structurally present, numerically zero
        }
      };
    }
  };

  // ---------------------------------------------------------------------------
  // Randomized mechanisms.
  // ---------------------------------------------------------------------------
  struct RandomReaction
  {
    std::vector<int> reactants;  // species indices, duplicates = order
    std::vector<int> products;
    double k;
  };

  struct RandomMechanism
  {
    int n_species;
    std::vector<RandomReaction> reactions;
    std::vector<int> constrained;  // QSSA species indices
  };

  std::string SpeciesName(int index)
  {
    return "S" + std::to_string(index);
  }

  // Net stoichiometry of species `s` in reaction `r`.
  int NetStoich(const RandomReaction& r, int s)
  {
    int net = 0;
    for (int p : r.products)
      net += (p == s);
    for (int q : r.reactants)
      net -= (q == s);
    return net;
  }

  RandomMechanism GenerateMechanism(unsigned seed)
  {
    std::mt19937 rng(seed);
    RandomMechanism mech;
    mech.n_species = 4 + static_cast<int>(rng() % 7);  // 4..10
    const int n_reactions = mech.n_species + static_cast<int>(rng() % mech.n_species);
    std::uniform_real_distribution<double> log_rate(-2.0, 2.0);
    auto random_species = [&]() { return static_cast<int>(rng() % mech.n_species); };

    for (int j = 0; j < n_reactions; ++j)
    {
      RandomReaction r;
      const int n_react = 1 + static_cast<int>(rng() % 2);
      const int n_prod = 1 + static_cast<int>(rng() % 2);
      for (int i = 0; i < n_react; ++i)
        r.reactants.push_back(random_species());
      for (int i = 0; i < n_prod; ++i)
        r.products.push_back(random_species());
      r.k = std::pow(10.0, log_rate(rng));
      mech.reactions.push_back(r);
    }

    // Choose 1-2 constrained species; forbid them from appearing twice as a
    // reactant in any reaction (keeps their constraint row linear in
    // themselves), then guarantee a strong linear loss so the row is index 1
    // and genuinely fast.
    const int n_constrained = 1 + static_cast<int>(rng() % 2);
    std::set<int> chosen;
    while (static_cast<int>(chosen.size()) < n_constrained)
      chosen.insert(random_species());
    mech.constrained.assign(chosen.begin(), chosen.end());

    for (auto& r : mech.reactions)
    {
      for (int c : mech.constrained)
      {
        int count = 0;
        for (int q : r.reactants)
          count += (q == c);
        while (count > 1)
        {
          for (auto& q : r.reactants)
          {
            if (q == c)
            {
              q = (c + 1) % mech.n_species;
              break;
            }
          }
          --count;
        }
      }
    }
    for (int c : mech.constrained)
    {
      RandomReaction loss;
      loss.reactants = { c };
      int sink = random_species();
      if (sink == c)
        sink = (c + 1) % mech.n_species;
      loss.products = { sink };
      loss.k = 1.0e4;  // dominant fast loss
      mech.reactions.push_back(loss);
    }
    return mech;
  }

  /// Generic QSSA constraint over the reaction table: G_c = net production of
  /// each constrained species; Jacobian assembled from mass-action derivatives.
  class RandomQssaModel
  {
   public:
    explicit RandomQssaModel(RandomMechanism mech)
        : mech_(std::move(mech))
    {
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      std::set<std::string> names;
      for (int c : mech_.constrained)
        names.insert(SpeciesName(c));
      return names;
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      std::set<std::string> names;
      for (int s = 0; s < mech_.n_species; ++s)
        names.insert(SpeciesName(s));
      return names;
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      std::set<std::pair<std::size_t, std::size_t>> elements;
      for (int c : mech_.constrained)
      {
        const auto row = s.at(SpeciesName(c));
        for (const auto& r : mech_.reactions)
        {
          if (NetStoich(r, c) == 0)
            continue;
          for (int q : r.reactants)
            elements.insert({ row, s.at(SpeciesName(q)) });
        }
        elements.insert({ row, row });
      }
      return elements;
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
      auto mech = mech_;
      std::vector<std::size_t> index_of(mech.n_species);
      for (int i = 0; i < mech.n_species; ++i)
        index_of[i] = s.at(SpeciesName(i));
      return [mech, index_of](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, DenseMatrixPolicy& f)
      {
        for (std::size_t cell = 0; cell < y.NumRows(); ++cell)
        {
          for (int c : mech.constrained)
          {
            double net = 0.0;
            for (const auto& r : mech.reactions)
            {
              const int stoich = NetStoich(r, c);
              if (stoich == 0)
                continue;
              double rate = r.k;
              for (int q : r.reactants)
                rate *= y[cell][index_of[q]];
              net += stoich * rate;
            }
            f[cell][index_of[c]] = net;
          }
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      auto mech = mech_;
      std::vector<std::size_t> index_of(mech.n_species);
      for (int i = 0; i < mech.n_species; ++i)
        index_of[i] = s.at(SpeciesName(i));
      return [mech, index_of](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
      {
        for (std::size_t block = 0; block < jac.NumberOfBlocks(); ++block)
        {
          for (int c : mech.constrained)
          {
            const auto row = index_of[c];
            for (const auto& r : mech.reactions)
            {
              const int stoich = NetStoich(r, c);
              if (stoich == 0)
                continue;
              for (std::size_t which = 0; which < r.reactants.size(); ++which)
              {
                double derivative = r.k;
                for (std::size_t other = 0; other < r.reactants.size(); ++other)
                {
                  if (other != which)
                    derivative *= y[block][index_of[r.reactants[other]]];
                }
                jac[block][row][index_of[r.reactants[which]]] -= stoich * derivative;
              }
            }
          }
        }
      };
    }

    /// Relative constraint satisfaction of `values`: |G_c| / (loss scale).
    double MaxRelativeResidual(const std::vector<double>& values) const
    {
      double worst = 0.0;
      for (int c : mech_.constrained)
      {
        double net = 0.0;
        double gross = 1.0e-300;
        for (const auto& r : mech_.reactions)
        {
          const int stoich = NetStoich(r, c);
          if (stoich == 0)
            continue;
          double rate = r.k;
          for (int q : r.reactants)
            rate *= values[q];
          net += stoich * rate;
          gross += std::abs(stoich) * rate;
        }
        // Below the tolerance-implied absolute floor (|net| ~ L * atol with
        // L ~ 1e4 and atol 1e-14) the relative measure is meaningless: species
        // driven to denormal scales satisfy the constraint absolutely.
        worst = std::max(worst, std::abs(net) / std::max(gross, 1.0e-9));
      }
      return worst;
    }

   private:
    RandomMechanism mech_;
  };
}  // namespace

TEST(StructuralDae, Index2FormulationFailsCleanly)
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
                    .AddExternalModel(NotIndexOneModel())
                    .SetReorderState(false)
                    .Build();
  auto state = solver.GetState(1);
  state.SetRelativeTolerance(1.0e-6);
  state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
  const double x0 = 1.0, z0 = 0.25;
  state.variables_[0][state.variable_map_.at("X")] = x0;
  state.variables_[0][state.variable_map_.at("Z")] = z0;
  state.conditions_[0].temperature_ = 298.0;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);

  auto result = solver.Solve(1.0, state);

  // Named failure, no silent nonsense, and the caller's state rolled back.
  EXPECT_NE(result.state_, SolverState::Converged);
  EXPECT_EQ(state.variables_[0][state.variable_map_.at("X")], x0);
  EXPECT_EQ(state.variables_[0][state.variable_map_.at("Z")], z0);
}

TEST(StructuralDae, RandomizedMechanismSweep)
{
  constexpr int kMechanisms = 200;
  int failures = 0;
  for (unsigned seed = 0; seed < kMechanisms; ++seed)
  {
    auto mech = GenerateMechanism(seed);
    RandomQssaModel model(mech);

    std::vector<PhaseSpecies> phase_species;
    for (int i = 0; i < mech.n_species; ++i)
      phase_species.push_back(Species(SpeciesName(i)));
    Phase gas_phase{ "gas", phase_species };

    std::vector<Process> processes;
    for (const auto& r : mech.reactions)
    {
      std::vector<Species> reactants;
      for (int q : r.reactants)
        reactants.push_back(Species(SpeciesName(q)));
      std::vector<StoichSpecies> products;
      for (int p : r.products)
        products.push_back({ Species(SpeciesName(p)), 1.0 });
      processes.push_back(ChemicalReactionBuilder()
                              .SetReactants(reactants)
                              .SetProducts(products)
                              .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = r.k })
                              .SetPhase(gas_phase)
                              .Build());
    }

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(std::move(options))
                      .SetSystem(System(gas_phase))
                      .SetReactions(processes)
                      .AddExternalModel(model)
                      .SetReorderState(false)
                      .Build();
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(1.0e-6);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-14));
    std::mt19937 ic_rng(seed ^ 0x9e3779b9u);
    std::uniform_real_distribution<double> ic(0.1, 1.0);
    for (int i = 0; i < mech.n_species; ++i)
      state.variables_[0][state.variable_map_.at(SpeciesName(i))] = ic(ic_rng);
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    bool converged = true;
    double advanced = 0.0;
    int guard = 0;
    while (advanced < 1.0 && guard++ < 10000)
    {
      auto result = solver.Solve(1.0 - advanced, state);
      if (result.state_ != SolverState::Converged && result.state_ != SolverState::ConvergenceExceededMaxSteps)
      {
        converged = false;
        break;
      }
      if (result.stats_.final_time_ <= 0.0)
      {
        converged = false;
        break;
      }
      advanced += result.stats_.final_time_;
    }

    bool healthy = converged;
    if (converged)
    {
      std::vector<double> values(mech.n_species);
      for (int i = 0; i < mech.n_species; ++i)
      {
        values[i] = state.variables_[0][state.variable_map_.at(SpeciesName(i))];
        // The post-clamp algebraic reprojection can land a hair below zero
        // (observed: -5e-43); anything beyond roundoff is a real failure.
        healthy = healthy && std::isfinite(values[i]) && values[i] >= -1.0e-12;
      }
      healthy = healthy && model.MaxRelativeResidual(values) < 1.0e-6;
    }
    if (!healthy)
    {
      ++failures;
      double residual = -1.0;
      double min_value = 0.0;
      if (converged)
      {
        std::vector<double> values(mech.n_species);
        min_value = 1e300;
        for (int i = 0; i < mech.n_species; ++i)
        {
          values[i] = state.variables_[0][state.variable_map_.at(SpeciesName(i))];
          min_value = std::min(min_value, values[i]);
        }
        residual = model.MaxRelativeResidual(values);
      }
      std::cout << "randomized mechanism failure at seed " << seed << " (converged=" << converged
                << ", residual=" << residual << ", min_value=" << min_value << ")\n";
    }
  }
  EXPECT_EQ(failures, 0) << failures << " of " << kMechanisms << " randomized mechanisms failed";
}
