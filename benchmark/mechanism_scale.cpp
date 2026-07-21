// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Mechanism scale-up along the production axis: grid cells. In a host
// atmosphere model (CAM-chem / MUSICA) the mechanism is fixed and the state
// carries many columns, so the question that decides production use is how
// the QSSA-DAE's cost scales with the cell count relative to the full ODE.
// This benchmark integrates the 16-species tropospheric mechanism at
// j_scale = 1 for 1..1000 identical cells and reports accepted steps and
// wall-clock per cell for both formulations (interleaved repetitions).
//
// Scale-up to CB05/TS1-class mechanisms (10-30 constrained radicals) requires
// mechanism-import infrastructure and is recorded as follow-on work in the
// plan; the cell axis is measured here because it is the one a 3-D host
// exercises hardest.
#include "tropospheric_qssa_constraint.hpp"
#include "tropospheric_system.hpp"

#include <micm/CPU.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  constexpr double T_BOX = 298.15;
  constexpr double P_BOX = 101325.0;
  constexpr double N_M = 2.5e19;
  constexpr double N_N2 = 0.78 * N_M;
  constexpr double N_O2 = 0.21 * N_M;
  constexpr double N_H2O = 0.01 * N_M;
  constexpr double N_O3_0 = 7.5e11;
  constexpr double N_NO_0 = 2.5e10;
  constexpr double N_NO2_0 = 2.5e10;
  constexpr double N_CO_0 = 2.5e12;
  constexpr double N_CH4_0 = 4.4e13;
  constexpr double N_CO2_0 = 1.0e16;

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  template<typename Solver>
  void InitState(Solver& solver, auto& state, double rtol, const tropospheric::ConsistentRadicals& rad0)
  {
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
    const std::size_t n_cells = state.number_of_grid_cells_;
    state.SetCustomRateParameter("p1", std::vector<double>(n_cells, tropospheric::JNO2_DEFAULT));
    state.SetCustomRateParameter("p2", std::vector<double>(n_cells, tropospheric::JO1D_DEFAULT));
    state.SetCustomRateParameter("p3", std::vector<double>(n_cells, tropospheric::JO3P_DEFAULT));
    auto& m = state.variable_map_;
    for (std::size_t cell = 0; cell < state.number_of_grid_cells_; ++cell)
    {
      auto set = [&](const std::string& name, double value) { state.variables_[cell][m.at(name)] = value; };
      set("O3", N_O3_0);
      set("NO", N_NO_0);
      set("NO2", N_NO2_0);
      set("CO", N_CO_0);
      set("CH4", N_CH4_0);
      set("CO2", N_CO2_0);
      set("H2O2", 0.0);
      set("HNO3", 0.0);
      set("O2", N_O2);
      set("N2", N_N2);
      set("M", N_M);
      set("H2O", N_H2O);
      set("O1D", rad0.O1D);
      set("O", rad0.O);
      set("OH", rad0.OH);
      set("HO2", rad0.HO2);
      state.conditions_[cell].temperature_ = T_BOX;
      state.conditions_[cell].pressure_ = P_BOX;
      state.conditions_[cell].air_density_ = N_M;
    }
    solver.UpdateStateParameters(state);
  }

  struct ScaleResult
  {
    std::uint64_t accepted = 0;
    double wall_us = 0.0;
    bool converged = true;
  };

  template<typename Solver>
  ScaleResult RunOnce(Solver& solver, std::size_t cells, double rtol, const std::vector<double>& output_times,
                      const tropospheric::ConsistentRadicals& rad0, bool timed_only)
  {
    auto state = solver.GetState(cells);
    InitState(solver, state, rtol, rad0);
    ScaleResult out;
    const auto t0 = std::chrono::steady_clock::now();
    double current = 0.0;
    for (double t_out : output_times)
    {
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
        {
          out.converged = false;
          return out;
        }
        if (result.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          return out;
        }
        if (!timed_only)
          out.accepted += result.stats_.accepted_;
        done += result.stats_.final_time_;
      }
      current = t_out;
    }
    const auto t1 = std::chrono::steady_clock::now();
    out.wall_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    return out;
  }
}  // namespace

int main()
{
  const double rtol = 1e-6;
  const int reps = 5;
  const auto output_times = OutputTimes(1e-3, 1e5, 10);
  const std::vector<std::size_t> cell_counts = { 1, 10, 100, 1000 };

  auto sys = tropospheric::MakeSystem();
  auto rates = tropospheric::RadicalRates::At(
      T_BOX, P_BOX, tropospheric::JNO2_DEFAULT, tropospheric::JO1D_DEFAULT, tropospheric::JO3P_DEFAULT);
  const auto rad0 = tropospheric::ProjectRadicals(rates, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O);
  const tropospheric::RefConc ref{ .O1D = rad0.O1D, .O = rad0.O,   .OH = rad0.OH, .HO2 = rad0.HO2,
                                   .O3 = N_O3_0,    .NO = N_NO_0,  .NO2 = N_NO2_0, .CO = N_CO_0,
                                   .CH4 = N_CH4_0,  .N2 = N_N2,    .O2 = N_O2,     .M = N_M,
                                   .H2O = N_H2O };
  tropospheric::QssaRadicalConstraint constraint(rates, ref);

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto ode_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(sys.gas_phase))
                        .SetReactions(sys.processes)
                        .SetReorderState(false)
                        .Build();
  auto dae_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(sys.gas_phase))
                        .SetReactions(sys.processes)
                        .AddExternalModel(constraint)
                        .SetReorderState(false)
                        .Build();

  std::ofstream csv("mechanism_scale.csv");
  csv << "cells,method,accepted,wallclock_median_us,us_per_cell,converged\n";
  csv.precision(12);

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "Tropospheric cell-count scaling (j_scale = 1, rtol = " << rtol << ", identical cells)\n";
  for (const std::size_t cells : cell_counts)
  {
    auto ode_stats = RunOnce(ode_solver, cells, rtol, output_times, rad0, false);
    auto dae_stats = RunOnce(dae_solver, cells, rtol, output_times, rad0, false);
    std::vector<double> ode_samples, dae_samples;
    for (int rep = 0; rep < reps; ++rep)
    {
      ode_samples.push_back(RunOnce(ode_solver, cells, rtol, output_times, rad0, true).wall_us);
      dae_samples.push_back(RunOnce(dae_solver, cells, rtol, output_times, rad0, true).wall_us);
    }
    std::sort(ode_samples.begin(), ode_samples.end());
    std::sort(dae_samples.begin(), dae_samples.end());
    const double ode_wall = ode_samples[ode_samples.size() / 2];
    const double dae_wall = dae_samples[dae_samples.size() / 2];
    csv << cells << ",full_ode," << ode_stats.accepted << ',' << ode_wall << ',' << ode_wall / cells << ','
        << (ode_stats.converged ? 1 : 0) << '\n';
    csv << cells << ",dae_qssa," << dae_stats.accepted << ',' << dae_wall << ',' << dae_wall / cells << ','
        << (dae_stats.converged ? 1 : 0) << '\n';
    std::cout << "cells=" << cells << "  ODE steps=" << ode_stats.accepted << " us/cell=" << ode_wall / cells
              << "  DAE steps=" << dae_stats.accepted << " us/cell=" << dae_wall / cells
              << "  ratio=" << dae_wall / ode_wall << (ode_stats.converged && dae_stats.converged ? "" : "  [FAILED]")
              << "\n";
  }
  std::cout << "wrote mechanism_scale.csv\n";
  return 0;
}
