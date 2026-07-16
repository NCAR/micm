// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Benchmark: full-ODE vs QSSA-DAE Robertson across a stiffness (k2) sweep.
// Emits deterministic solver-stat counters (primary evidence), median
// wall-clock (secondary, hardware-dependent), and post-transient accuracy of
// the DAE vs a tight-tolerance full-ODE reference. Writes one CSV.
#include "robertson_qssa_constraint.hpp"
#include "robertson_system.hpp"

#include <micm/CPU.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
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

  // Aggregated solver counters plus the trajectory sampled at the output times.
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
    std::vector<double> a_at_output;  // one entry per output time
    std::vector<double> c_at_output;
  };

  // Logarithmically spaced output times in [t_first, t_last].
  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    double log_first = std::log10(t_first);
    double log_last = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, log_first + (log_last - log_first) * i / (n - 1)));
    return times;
  }

  // Run one integration to completion, accumulating solver stats across every
  // Solve() call and recording A and C at each output time. The wallclock_reps
  // passes are timed only (stats ignored) for a median wall-clock estimate.
  CaseResult
  RunCase(Method method, double k1, double k2, double k3, double rtol, const std::vector<double>& output_times, int wallclock_reps)
  {
    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    // Build in one unbroken chain (never store the intermediate builder in
    // `auto` — chained setters return a base SolverBuilder& and would slice).
    // Both branches yield the SAME solver type (AddExternalModel returns
    // SolverBuilder&, type-erased), so the ternary is well-typed.
    robertson::QssaConstraint constraint(k1, k2, k3);

    auto solver =
        (method == Method::DaeQssa)
            ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(sys.gas_phase))
                  .SetReactions(sys.processes)
                  .AddExternalModel(constraint)
                  .SetReorderState(false)
                  .Build()
            : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(sys.gas_phase))
                  .SetReactions(sys.processes)
                  .SetReorderState(false)
                  .Build();

    auto init_state = [&](auto& state)
    {
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, rtol * 1e-2));
      state.SetCustomRateParameter("r1", k1);
      state.SetCustomRateParameter("r2", k2);
      state.SetCustomRateParameter("r3", k3);
      auto map = state.variable_map_;
      state.variables_[0][map.at("A")] = 1.0;
      state.variables_[0][map.at("B")] = (method == Method::DaeQssa) ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
      state.variables_[0][map.at("C")] = 0.0;
      state.conditions_[0].temperature_ = 272.5;
      state.conditions_[0].pressure_ = 101253.3;
      state.conditions_[0].air_density_ = 1e6;
      solver.UpdateStateParameters(state);
    };

    CaseResult out;

    // ---- Stats + trajectory pass (untimed) ----
    {
      auto state = solver.GetState(1);
      init_state(state);
      auto map = state.variable_map_;
      double current = 0.0;
      for (double t_out : output_times)
      {
        double done = 0.0;
        double dt = t_out - current;
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
        out.a_at_output.push_back(state.variables_[0][map.at("A")]);
        out.c_at_output.push_back(state.variables_[0][map.at("C")]);
        current = t_out;
        if (!out.converged)
          break;
      }
    }

    // ---- Wall-clock passes (timed, stats ignored) ----
    std::vector<double> samples;
    for (int rep = 0; rep < wallclock_reps; ++rep)
    {
      auto state = solver.GetState(1);
      init_state(state);
      auto t0 = std::chrono::steady_clock::now();
      double current = 0.0;
      for (double t_out : output_times)
      {
        double done = 0.0;
        double dt = t_out - current;
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
      auto t1 = std::chrono::steady_clock::now();
      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    if (!samples.empty())
    {
      std::sort(samples.begin(), samples.end());
      out.wallclock_median_us = samples[samples.size() / 2];
    }
    return out;
  }

  // Max relative error of A and C of `cand` vs `ref`, over output times >= t_skip.
  double
  MaxRelErrorPostTransient(const CaseResult& cand, const CaseResult& ref, const std::vector<double>& output_times, double t_skip)
  {
    double worst = 0.0;
    for (std::size_t i = 0; i < output_times.size(); ++i)
    {
      if (output_times[i] < t_skip)
        continue;
      if (i >= cand.a_at_output.size() || i >= ref.a_at_output.size())
        break;
      auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
      worst = std::max(worst, rel(cand.a_at_output[i], ref.a_at_output[i]));
      worst = std::max(worst, rel(cand.c_at_output[i], ref.c_at_output[i]));
    }
    return worst;
  }
}  // namespace

int main()
{
  const double k1 = robertson::K1_DEFAULT;
  const double k3 = robertson::K3_DEFAULT;
  const double rtol = 1e-6;
  const double t_skip = 1.0;     // transient cutoff (C's QSSA transient outlasts B's by ~2 decades)
  const int wallclock_reps = 7;  // odd -> well-defined median
  const std::vector<double> k2_sweep = { 3e5, 3e6, 3e7, 3e8, 3e9 };
  const auto output_times = OutputTimes(1e-3, 1e6, 19);  // spans transient + long tail

  std::ofstream csv("robertson_dae_benchmark.csv");
  csv << "method,k2,rtol,number_of_steps,accepted,rejected,function_calls,"
         "jacobian_updates,decompositions,solves,converged,wallclock_median_us,max_rel_err\n";

  auto write_row = [&](const std::string& method, double k2, const CaseResult& r, double max_rel_err)
  {
    csv << method << ',' << k2 << ',' << rtol << ',' << r.number_of_steps << ',' << r.accepted << ',' << r.rejected << ','
        << r.function_calls << ',' << r.jacobian_updates << ',' << r.decompositions << ',' << r.solves << ','
        << (r.converged ? 1 : 0) << ',' << r.wallclock_median_us << ',' << max_rel_err << '\n';
  };

  std::cout << "Robertson: full-ODE vs QSSA-DAE across k2 (rtol=" << rtol << ")\n";
  for (double k2 : k2_sweep)
  {
    auto reference = RunCase(Method::FullOde, k1, k2, k3, 1e-10, output_times, 0);  // tight-tol accuracy reference
    auto ode = RunCase(Method::FullOde, k1, k2, k3, rtol, output_times, wallclock_reps);
    auto dae = RunCase(Method::DaeQssa, k1, k2, k3, rtol, output_times, wallclock_reps);

    double ode_err = MaxRelErrorPostTransient(ode, reference, output_times, t_skip);
    double dae_err = MaxRelErrorPostTransient(dae, reference, output_times, t_skip);

    write_row("full_ode", k2, ode, ode_err);
    write_row("dae_qssa", k2, dae, dae_err);

    std::cout << "k2=" << k2 << "  ODE steps=" << ode.number_of_steps << " (us=" << ode.wallclock_median_us << ")"
              << "  DAE steps=" << dae.number_of_steps << " (us=" << dae.wallclock_median_us << ")"
              << "  DAE max_rel_err=" << dae_err << "\n";
  }

  std::cout << "wrote robertson_dae_benchmark.csv\n";
  return 0;
}
