// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Schur-reduction ceiling measurement (phase 4 prototype). The equilibrium
// family is the case where the DAE loses on wall-clock because MICM keeps
// algebraic variables in the factored linear system: 3N rows with two-thirds
// algebraic. A production Schur reduction would factor only the differential
// block. This prototype measures the CEILING of that optimization by solving
// the exactly reduced problem: with B_i = Keq*A_i and A_i + B_i + C_i = Ctot
// eliminated analytically, the slow dynamics collapse to N independent linear
// ODEs D_i' = -a D_i with D = Ctot - C and a = ks*Keq/(1+Keq).
//
//   full_ode    - 3N stiff ODE resolving the S-fast equilibration
//   dae         - 3N with equilibrium + conservation constraints (A, B algebraic)
//   reduced_ode - N-dimensional exact reduction (the Schur ceiling)
//
// The reduced model starts on the slow manifold (it has no fast transient),
// so wall-clock is compared on the shared slow horizon; the ceiling ratio
// reduced/full quantifies the maximum achievable gain for the core
// implementation designed in
// docs/superpowers/notes/2026-07-20-schur-reduction-design.md.
#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace micm;

namespace
{
  constexpr double Keq = 1.0;
  constexpr double Ctot = 1.0;
  constexpr double ks = 1.0;
  constexpr double S = 1.0e6;
  constexpr double T_FINAL = 20.0;

  struct Result
  {
    std::uint64_t steps = 0;
    std::uint64_t rejected = 0;
    double wall_us = 0.0;
    bool converged = true;
  };

  template<typename Solver>
  Result RunLoop(Solver& solver, auto init, int reps)
  {
    // One state reused across repetitions: initialization resets the
    // variables, and any solver-side caches (e.g. the Schur symbolic
    // structures) warm up in the stats pass — the production pattern, where a
    // host model owns long-lived states.
    auto state = solver.GetState(1);
    Result out;
    {
      init(state);
      solver.UpdateStateParameters(state);
      double done = 0.0;
      while (done < T_FINAL)
      {
        auto r = solver.Solve(T_FINAL - done, state);
        if ((r.state_ != SolverState::Converged && r.state_ != SolverState::ConvergenceExceededMaxSteps) ||
            r.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          return out;
        }
        out.steps += r.stats_.accepted_;
        out.rejected += r.stats_.rejected_;
        done += r.stats_.final_time_;
      }
    }
    std::vector<double> samples;
    for (int rep = 0; rep < reps; ++rep)
    {
      init(state);
      solver.UpdateStateParameters(state);
      const auto t0 = std::chrono::steady_clock::now();
      double done = 0.0;
      while (done < T_FINAL)
      {
        auto r = solver.Solve(T_FINAL - done, state);
        if ((r.state_ != SolverState::Converged && r.state_ != SolverState::ConvergenceExceededMaxSteps) ||
            r.stats_.final_time_ <= 0.0)
          break;
        done += r.stats_.final_time_;
      }
      const auto t1 = std::chrono::steady_clock::now();
      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    std::sort(samples.begin(), samples.end());
    out.wall_us = samples[samples.size() / 2];
    return out;
  }

  auto MakeCommonInit(int N, bool dae)
  {
    return [N, dae](auto& state)
    {
      state.SetRelativeTolerance(1.0e-6);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
      auto& vm = state.variable_map_;
      for (int i = 0; i < N; ++i)
      {
        // Off-equilibrium start: all mass in A. The DAE projects onto the
        // manifold at initialization; the ODE resolves the S-fast transient.
        state.variables_[0][vm.at("A" + std::to_string(i))] = dae ? Ctot / (1.0 + Keq) : Ctot;
        state.variables_[0][vm.at("B" + std::to_string(i))] = dae ? Ctot * Keq / (1.0 + Keq) : 0.0;
        state.variables_[0][vm.at("C" + std::to_string(i))] = 0.0;
      }
      state.conditions_[0].temperature_ = 298.0;
      state.conditions_[0].pressure_ = 101325.0;
    };
  }
}  // namespace

int main()
{
  const std::vector<int> sizes = { 1, 4, 16, 64, 256 };
  const int reps = 7;

  std::ofstream csv("schur_ceiling.csv");
  csv << "N,method,accepted,wallclock_median_us\n";
  csv.precision(12);

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "Schur ceiling: full 3N ODE vs 3N DAE vs exact N-dim reduction (S=" << S << ")\n";
  for (const int N : sizes)
  {
    // ---- full ODE and DAE share the 3N system ----
    std::vector<PhaseSpecies> species;
    for (int i = 0; i < N; ++i)
    {
      species.push_back(Species("A" + std::to_string(i)));
      species.push_back(Species("B" + std::to_string(i)));
      species.push_back(Species("C" + std::to_string(i)));
    }
    Phase gas{ "gas", species };

    std::vector<Process> ode_processes;
    std::vector<Process> slow_processes;
    for (int i = 0; i < N; ++i)
    {
      auto A = Species("A" + std::to_string(i));
      auto B = Species("B" + std::to_string(i));
      auto C = Species("C" + std::to_string(i));
      ode_processes.push_back(ChemicalReactionBuilder()
                                  .SetReactants({ A })
                                  .SetProducts({ { B, 1 } })
                                  .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = S * Keq })
                                  .SetPhase(gas)
                                  .Build());
      ode_processes.push_back(ChemicalReactionBuilder()
                                  .SetReactants({ B })
                                  .SetProducts({ { A, 1 } })
                                  .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = S })
                                  .SetPhase(gas)
                                  .Build());
      Process slow = ChemicalReactionBuilder()
                         .SetReactants({ B })
                         .SetProducts({ { C, 1 } })
                         .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = ks })
                         .SetPhase(gas)
                         .Build();
      ode_processes.push_back(slow);
      slow_processes.push_back(slow);
    }

    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto ode_solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                          .SetSystem(System(gas))
                          .SetReactions(ode_processes)
                          .SetReorderState(false)
                          .Build();

    std::vector<Constraint> constraints;
    for (int i = 0; i < N; ++i)
    {
      auto A = Species("A" + std::to_string(i));
      auto B = Species("B" + std::to_string(i));
      auto C = Species("C" + std::to_string(i));
      constraints.push_back(EquilibriumConstraint(
          "eq" + std::to_string(i),
          B,
          std::vector<StoichSpecies>{ { A, 1.0 } },
          std::vector<StoichSpecies>{ { B, 1.0 } },
          VantHoffParam{ .K_HLC_ref_ = Keq, .delta_H_ = 0.0 }));
      constraints.push_back(
          LinearConstraint("mass" + std::to_string(i), A, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, Ctot));
    }
    auto dae_solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                          .SetSystem(System(gas))
                          .SetReactions(slow_processes)
                          .SetConstraints(std::move(constraints))
                          .SetReorderState(false)
                          .Build();

    // ---- exact N-dimensional reduction: D_i' = -a D_i ----
    const double a = ks * Keq / (1.0 + Keq);
    std::vector<PhaseSpecies> reduced_species;
    for (int i = 0; i < N; ++i)
      reduced_species.push_back(Species("D" + std::to_string(i)));
    Phase reduced_gas{ "gas", reduced_species };
    std::vector<Process> reduced_processes;
    for (int i = 0; i < N; ++i)
    {
      reduced_processes.push_back(ChemicalReactionBuilder()
                                      .SetReactants({ Species("D" + std::to_string(i)) })
                                      .SetProducts({})
                                      .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = a })
                                      .SetPhase(reduced_gas)
                                      .Build());
    }
    auto reduced_solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                              .SetSystem(System(reduced_gas))
                              .SetReactions(reduced_processes)
                              .SetReorderState(false)
                              .Build();
    auto reduced_init = [N](auto& state)
    {
      state.SetRelativeTolerance(1.0e-6);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
      auto& vm = state.variable_map_;
      for (int i = 0; i < N; ++i)
        state.variables_[0][vm.at("D" + std::to_string(i))] = Ctot;
      state.conditions_[0].temperature_ = 298.0;
      state.conditions_[0].pressure_ = 101325.0;
    };

    std::vector<Constraint> schur_constraints;
    for (int i = 0; i < N; ++i)
    {
      auto A = Species("A" + std::to_string(i));
      auto B = Species("B" + std::to_string(i));
      auto C = Species("C" + std::to_string(i));
      schur_constraints.push_back(EquilibriumConstraint(
          "eq" + std::to_string(i),
          B,
          std::vector<StoichSpecies>{ { A, 1.0 } },
          std::vector<StoichSpecies>{ { B, 1.0 } },
          VantHoffParam{ .K_HLC_ref_ = Keq, .delta_H_ = 0.0 }));
      schur_constraints.push_back(
          LinearConstraint("mass" + std::to_string(i), A, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, Ctot));
    }
    auto schur_options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    schur_options.schur_reduction_ = true;
    auto schur_solver = CpuSolverBuilder<RosenbrockSolverParameters>(schur_options)
                            .SetSystem(System(gas))
                            .SetReactions(slow_processes)
                            .SetConstraints(std::move(schur_constraints))
                            .SetReorderState(false)
                            .Build();

    auto ode = RunLoop(ode_solver, MakeCommonInit(N, false), reps);
    auto dae = RunLoop(dae_solver, MakeCommonInit(N, true), reps);
    auto dae_schur = RunLoop(schur_solver, MakeCommonInit(N, true), reps);
    auto reduced = RunLoop(reduced_solver, reduced_init, reps);

    csv << N << ",full_ode," << ode.steps << ',' << ode.wall_us << '\n';
    csv << N << ",dae," << dae.steps << ',' << dae.wall_us << '\n';
    csv << N << ",dae_schur," << dae_schur.steps << ',' << dae_schur.wall_us << '\n';
    csv << N << ",reduced_ode," << reduced.steps << ',' << reduced.wall_us << '\n';
    std::cout << "N=" << N << "  ODE us=" << ode.wall_us << " (acc=" << ode.steps << " rej=" << ode.rejected
              << ")  DAE us=" << dae.wall_us << " (acc=" << dae.steps << " rej=" << dae.rejected
              << ")  DAE+Schur us=" << dae_schur.wall_us << "  reduced us=" << reduced.wall_us << " (acc="
              << reduced.steps << ")  DAE/ODE=" << dae.wall_us / ode.wall_us
              << "  Schur/ODE=" << dae_schur.wall_us / ode.wall_us
              << "  ceiling reduced/ODE=" << reduced.wall_us / ode.wall_us
              << (ode.converged && dae.converged && dae_schur.converged && reduced.converged ? "" : "  [FAILED]")
              << "\n";
  }
  std::cout << "wrote schur_ceiling.csv\n";
  return 0;
}
