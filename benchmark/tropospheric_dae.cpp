// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Benchmark: full-ODE vs QSSA-DAE for the tropospheric O3-NOx-HOx mechanism,
// swept over a fixed (steady) photolysis intensity. Each run holds J constant;
// the sweep scales jNO2/jO1D/jO3P by a single photolysis-scale factor (j_scale). This is the
// higher-dimensional successor to robertson_dae: FOUR fast radicals (O1D, O,
// OH, HO2) are removed by QSSA in the DAE branch instead of one species.
//
// The full ODE starts all radicals at zero and must resolve the radical
// transients (O1D ~ ns up to HO2 ~ seconds); the DAE starts on the QSSA
// manifold and skips them. We report deterministic solver-stat counters
// (primary evidence), median wall-clock (secondary), and post-transient
// accuracy of slow observables (O3, HNO3) vs a tight-tolerance ODE reference.
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
  enum class Method
  {
    FullOde,
    DaeQssa
  };

  // Box conditions and initial slow-species number densities (molec/cm^3),
  // representative of the polluted continental boundary layer.
  constexpr double T_BOX = 298.15;     // K
  constexpr double P_BOX = 101325.0;   // Pa
  constexpr double N_M = 2.5e19;       // total air number density
  constexpr double N_N2 = 0.78 * N_M;
  constexpr double N_O2 = 0.21 * N_M;
  constexpr double N_H2O = 0.01 * N_M;  // ~1% H2O
  constexpr double N_O3_0 = 7.5e11;     // ~30 ppb
  constexpr double N_NO_0 = 2.5e10;     // ~1 ppb
  constexpr double N_NO2_0 = 2.5e10;    // ~1 ppb
  constexpr double N_CO_0 = 2.5e12;     // ~100 ppb
  constexpr double N_CH4_0 = 4.4e13;    // ~1.8 ppm
  constexpr double N_CO2_0 = 1.0e16;    // bath (product sink target)

  struct CaseResult
  {
    std::uint64_t number_of_steps = 0;
    std::uint64_t accepted = 0;
    std::uint64_t rejected = 0;
    std::uint64_t function_calls = 0;
    std::uint64_t jacobian_updates = 0;
    std::uint64_t decompositions = 0;
    std::uint64_t solves = 0;
    bool converged = true;
    double wallclock_median_us = 0.0;
    std::vector<double> o3_at_output;
    std::vector<double> hno3_at_output;
  };

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  CaseResult RunCase(Method method, double j_scale, double rtol, const std::vector<double>& output_times, int wallclock_reps)
  {
    auto sys = tropospheric::MakeSystem();
    const double j1 = j_scale * tropospheric::JNO2_DEFAULT;
    const double j2 = j_scale * tropospheric::JO1D_DEFAULT;
    const double j3 = j_scale * tropospheric::JO3P_DEFAULT;
    auto rates = tropospheric::RadicalRates::At(T_BOX, P_BOX, j1, j2, j3);
    // Reference radicals (for per-row residual scaling) from the QSSA manifold
    // at the initial slow-species state.
    const auto ref_rad = tropospheric::ProjectRadicals(rates, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O);
    const tropospheric::RefConc ref{ .O1D = ref_rad.O1D, .O = ref_rad.O,   .OH = ref_rad.OH, .HO2 = ref_rad.HO2,
                                     .O3 = N_O3_0,       .NO = N_NO_0,     .NO2 = N_NO2_0,   .CO = N_CO_0,
                                     .CH4 = N_CH4_0,     .N2 = N_N2,       .O2 = N_O2,       .M = N_M,
                                     .H2O = N_H2O };
    tropospheric::QssaRadicalConstraint constraint(rates, ref);

    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver =
        (method == Method::DaeQssa)
            ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                  .SetReactions(sys.processes)
                  .AddExternalModel(constraint)
                  .SetReorderState(false)
                  .Build()
            : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                  .SetReactions(sys.processes)
                  .SetReorderState(false)
                  .Build();

    auto init_state = [&](auto& state)
    {
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));  // molec/cm^3 floor
      state.SetCustomRateParameter("p1", j1);
      state.SetCustomRateParameter("p2", j2);
      state.SetCustomRateParameter("p3", j3);
      auto& m = state.variable_map_;
      auto set = [&](const std::string& name, double val) { state.variables_[0][m.at(name)] = val; };
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
      // Both methods start on the QSSA manifold so they solve the SAME IVP:
      // this isolates the QSSA approximation error (and the solver-efficiency
      // comparison) from any radical spin-up IC difference. (A cold-start ODE
      // with radicals=0 solves a different IVP — its extra HNO3 deficit during
      // spin-up is permanent — so it is not a fair accuracy reference.)
      {
        auto c = tropospheric::ProjectRadicals(rates, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O);
        set("O1D", c.O1D);
        set("O", c.O);
        set("OH", c.OH);
        set("HO2", c.HO2);
      }
      state.conditions_[0].temperature_ = T_BOX;
      state.conditions_[0].pressure_ = P_BOX;
      state.conditions_[0].air_density_ = N_M;
      solver.UpdateStateParameters(state);
    };

    CaseResult out;

    // ---- Stats + trajectory pass (untimed) ----
    {
      auto state = solver.GetState(1);
      init_state(state);
      auto& m = state.variable_map_;
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
            break;
          }
          out.number_of_steps += result.stats_.number_of_steps_;
          out.accepted += result.stats_.accepted_;
          out.rejected += result.stats_.rejected_;
          out.function_calls += result.stats_.function_calls_;
          out.jacobian_updates += result.stats_.jacobian_updates_;
          out.decompositions += result.stats_.decompositions_;
          out.solves += result.stats_.solves_;
          if (result.stats_.final_time_ <= 0.0)
          {
            out.converged = false;
            break;
          }
          done += result.stats_.final_time_;
        }
        out.o3_at_output.push_back(state.variables_[0][m.at("O3")]);
        out.hno3_at_output.push_back(state.variables_[0][m.at("HNO3")]);
        current = t_out;
        if (!out.converged)
          break;
      }
    }

    // ---- Wall-clock passes (timed) ----
    std::vector<double> samples;
    for (int rep = 0; rep < wallclock_reps; ++rep)
    {
      auto state = solver.GetState(1);
      init_state(state);
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
            break;
          if (result.stats_.final_time_ <= 0.0)
            break;
          done += result.stats_.final_time_;
        }
        current = t_out;
      }
      const auto t1 = std::chrono::steady_clock::now();
      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    if (!samples.empty())
    {
      std::sort(samples.begin(), samples.end());
      out.wallclock_median_us = samples[samples.size() / 2];
    }
    return out;
  }

  double MaxRelErrorPostTransient(const CaseResult& cand, const CaseResult& ref, const std::vector<double>& output_times, double t_skip)
  {
    double worst = 0.0;
    for (std::size_t i = 0; i < output_times.size(); ++i)
    {
      if (output_times[i] < t_skip)
        continue;
      if (i >= cand.o3_at_output.size() || i >= ref.o3_at_output.size())
        break;
      auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
      worst = std::max(worst, rel(cand.o3_at_output[i], ref.o3_at_output[i]));
      worst = std::max(worst, rel(cand.hno3_at_output[i], ref.hno3_at_output[i]));
    }
    return worst;
  }
}  // namespace

int main()
{
  const double rtol = 1e-6;
  const double t_skip = 10.0;     // skip the radical transient (slowest radical HO2 ~ 5 s)
  const int wallclock_reps = 7;
  const std::vector<double> j_scale_sweep = { 0.1, 0.25, 0.5, 1.0, 2.0 };
  const auto output_times = OutputTimes(1e-3, 1e5, 25);  // ns transient up to ~1 day relaxation

  std::ofstream csv("tropospheric_dae_benchmark.csv");
  csv << "method,j_scale,rtol,number_of_steps,accepted,rejected,function_calls,"
         "jacobian_updates,decompositions,solves,converged,wallclock_median_us,max_rel_err\n";

  auto write_row = [&](const std::string& method, double j_scale, const CaseResult& r, double err)
  {
    csv << method << ',' << j_scale << ',' << rtol << ',' << r.number_of_steps << ',' << r.accepted << ',' << r.rejected << ','
        << r.function_calls << ',' << r.jacobian_updates << ',' << r.decompositions << ',' << r.solves << ','
        << (r.converged ? 1 : 0) << ',' << r.wallclock_median_us << ',' << err << '\n';
  };

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "Tropospheric O3-NOx-HOx: full-ODE vs QSSA-DAE across fixed photolysis scale j_scale (rtol=" << rtol << ")\n";
  for (double j_scale : j_scale_sweep)
  {
    auto reference = RunCase(Method::FullOde, j_scale, 1e-10, output_times, 0);
    auto ode = RunCase(Method::FullOde, j_scale, rtol, output_times, wallclock_reps);
    auto dae = RunCase(Method::DaeQssa, j_scale, rtol, output_times, wallclock_reps);

    const double ode_err = MaxRelErrorPostTransient(ode, reference, output_times, t_skip);
    const double dae_err = MaxRelErrorPostTransient(dae, reference, output_times, t_skip);

    write_row("full_ode", j_scale, ode, ode_err);
    write_row("dae_qssa", j_scale, dae, dae_err);

    std::cout << "j_scale=" << j_scale << "  ODE acc=" << ode.accepted << " (us=" << ode.wallclock_median_us << ")"
              << "  DAE acc=" << dae.accepted << " (us=" << dae.wallclock_median_us << ")"
              << "  DAE max_rel_err=" << dae_err << (dae.converged ? "" : "  [DAE FAILED]") << "\n";
  }

  std::cout << "wrote tropospheric_dae_benchmark.csv\n";
  return 0;
}
