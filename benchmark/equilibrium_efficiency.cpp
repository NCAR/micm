// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Efficiency hunt: does an algebraic-constraint DAE beat the full ODE on
// WALL-CLOCK for the implicit Rosenbrock solver, and does any advantage grow
// with system size?
//
// N independent fast-equilibrium triples (i = 1..N):
//   A_i <-> B_i   fast reversible equilibrium (rate scale S),  B_i/A_i -> Keq
//   B_i -> C_i    slow product formation (rate ks)
//
// Full ODE: 3N reactions, all 3N species differential.
// DAE: N slow reactions (B_i->C_i) + N equilibrium constraints (B_i = Keq*A_i)
//      + N conservation constraints (A_i+B_i+C_i = Ctot). A_i, B_i algebraic;
//      C_i differential. (Note: MICM keeps algebraic vars in the matrix, so the
//      linear-system dimension is 3N either way.)
//
// Off-equilibrium start (all mass in A_i) so the ODE must resolve the fast
// transient the DAE skips. We report accepted steps and median wall-clock vs N.
#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
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
  constexpr double S = 1.0e6;  // fast-equilibrium rate scale

  struct Result
  {
    std::uint64_t steps = 0;
    double wall_us = 0.0;
    bool converged = true;
  };

  std::vector<Species> MakeSpeciesList(int N)
  {
    std::vector<Species> sp;
    for (int i = 0; i < N; ++i)
    {
      sp.push_back(Species("A" + std::to_string(i)));
      sp.push_back(Species("B" + std::to_string(i)));
      sp.push_back(Species("C" + std::to_string(i)));
    }
    return sp;
  }

  template<typename Solver>
  Result RunLoop(Solver& solver, int N, int wallclock_reps)
  {
    auto setup = [&](auto& state)
    {
      auto& vm = state.variable_map_;
      for (int i = 0; i < N; ++i)
      {
        state.variables_[0][vm.at("A" + std::to_string(i))] = Ctot;  // off-equilibrium: all in A
        state.variables_[0][vm.at("B" + std::to_string(i))] = 0.0;
        state.variables_[0][vm.at("C" + std::to_string(i))] = 0.0;
      }
      state.conditions_[0].temperature_ = 298.0;
      state.conditions_[0].pressure_ = 101325.0;
      solver.UpdateStateParameters(state);
    };

    const double T = 20.0;
    Result out;
    {
      auto state = solver.GetState(1);
      state.SetRelativeTolerance(1.0e-6);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
      setup(state);
      double done = 0.0;
      while (done < T)
      {
        auto r = solver.Solve(T - done, state);
        if ((r.state_ != SolverState::Converged && r.state_ != SolverState::ConvergenceExceededMaxSteps) ||
            r.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          break;
        }
        out.steps += r.stats_.accepted_;
        done += r.stats_.final_time_;
      }
    }

    std::vector<double> samples;
    for (int rep = 0; rep < wallclock_reps; ++rep)
    {
      auto state = solver.GetState(1);
      state.SetRelativeTolerance(1.0e-6);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
      setup(state);
      auto t0 = std::chrono::steady_clock::now();
      double done = 0.0;
      while (done < T)
      {
        auto r = solver.Solve(T - done, state);
        if ((r.state_ != SolverState::Converged && r.state_ != SolverState::ConvergenceExceededMaxSteps) ||
            r.stats_.final_time_ <= 0.0)
          break;
        done += r.stats_.final_time_;
      }
      auto t1 = std::chrono::steady_clock::now();
      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    std::sort(samples.begin(), samples.end());
    out.wall_us = samples[samples.size() / 2];
    return out;
  }

  Result RunOde(int N, int reps)
  {
    auto sp = MakeSpeciesList(N);
    Phase gas{ "gas", std::vector<PhaseSpecies>(sp.begin(), sp.end()) };
    std::vector<Process> rxns;
    for (int i = 0; i < N; ++i)
    {
      auto A = Species("A" + std::to_string(i)), B = Species("B" + std::to_string(i)), C = Species("C" + std::to_string(i));
      rxns.push_back(ChemicalReactionBuilder().SetReactants({ A }).SetProducts({ { B, 1 } }).SetRateConstant(ArrheniusRateConstantParameters{ .A_ = S, .B_ = 0, .C_ = 0 }).SetPhase(gas).Build());
      rxns.push_back(ChemicalReactionBuilder().SetReactants({ B }).SetProducts({ { A, 1 } }).SetRateConstant(ArrheniusRateConstantParameters{ .A_ = S / Keq, .B_ = 0, .C_ = 0 }).SetPhase(gas).Build());
      rxns.push_back(ChemicalReactionBuilder().SetReactants({ B }).SetProducts({ { C, 1 } }).SetRateConstant(ArrheniusRateConstantParameters{ .A_ = ks, .B_ = 0, .C_ = 0 }).SetPhase(gas).Build());
    }
    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                      .SetSystem(System(gas))
                      .SetReactions(rxns)
                      .SetReorderState(false)
                      .Build();
    return RunLoop(solver, N, reps);
  }

  Result RunDae(int N, int reps)
  {
    auto sp = MakeSpeciesList(N);
    Phase gas{ "gas", std::vector<PhaseSpecies>(sp.begin(), sp.end()) };
    std::vector<Process> rxns;
    std::vector<Constraint> constraints;
    for (int i = 0; i < N; ++i)
    {
      auto A = Species("A" + std::to_string(i)), B = Species("B" + std::to_string(i)), C = Species("C" + std::to_string(i));
      rxns.push_back(ChemicalReactionBuilder().SetReactants({ B }).SetProducts({ { C, 1 } }).SetRateConstant(ArrheniusRateConstantParameters{ .A_ = ks, .B_ = 0, .C_ = 0 }).SetPhase(gas).Build());
      constraints.push_back(EquilibriumConstraint(
          "eq" + std::to_string(i), B, std::vector<StoichSpecies>{ { A, 1.0 } }, std::vector<StoichSpecies>{ { B, 1.0 } }, VantHoffParam{ .K_HLC_ref_ = Keq, .delta_H_ = 0.0 }));
      constraints.push_back(LinearConstraint(
          "mass" + std::to_string(i), A, { { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, Ctot));
    }
    auto options = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(options)
                      .SetSystem(System(gas))
                      .SetReactions(rxns)
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();
    return RunLoop(solver, N, reps);
  }
}  // namespace

int main()
{
  const int reps = 7;
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "N fast-equilibrium triples, off-equilibrium start, S=" << S << ". Wall-clock vs system size.\n\n";
  std::cout << "N      ODE_steps  ODE_us       DAE_steps  DAE_us       DAE/ODE_wall  conv\n";
  for (int N : { 1, 4, 16, 64, 256 })
  {
    auto ode = RunOde(N, reps);
    auto dae = RunDae(N, reps);
    std::cout << N << "      " << ode.steps << "        " << ode.wall_us << "    " << dae.steps << "        " << dae.wall_us
              << "    " << dae.wall_us / ode.wall_us << "      " << (ode.converged && dae.converged ? "ok" : "FAIL") << "\n";
  }
  return 0;
}
