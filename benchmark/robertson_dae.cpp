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
#include <sstream>
#include <string>
#include <vector>

#ifndef MICM_BENCHMARK_DATA_DIR
#define MICM_BENCHMARK_DATA_DIR "benchmark/data"
#endif

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

  // External reference (independent solver), A and C at the benchmark output times.
  struct ExternalReference
  {
    std::vector<double> times;
    std::vector<double> a;
    std::vector<double> c;
  };

  ExternalReference LoadExternalReference(const std::string& path)
  {
    ExternalReference ref;
    std::ifstream file(path);
    if (!file)
      return ref;
    std::string line;
    std::getline(file, line);  // header
    std::vector<std::string> names;
    {
      std::istringstream header(line);
      std::string field;
      while (std::getline(header, field, ','))
        names.push_back(field);
    }
    auto column = [&](const std::string& name) -> std::ptrdiff_t
    {
      for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i] == name)
          return static_cast<std::ptrdiff_t>(i);
      return -1;
    };
    const auto i_time = column("time"), i_a = column("A"), i_c = column("C");
    if (i_time < 0 || i_a < 0 || i_c < 0)
      return ref;
    while (std::getline(file, line))
    {
      std::istringstream row(line);
      std::string field;
      std::vector<double> values;
      while (std::getline(row, field, ','))
        values.push_back(std::stod(field));
      if (static_cast<std::ptrdiff_t>(values.size()) <= std::max({ i_time, i_a, i_c }))
        continue;
      ref.times.push_back(values[i_time]);
      ref.a.push_back(values[i_a]);
      ref.c.push_back(values[i_c]);
    }
    return ref;
  }

  double MaxRelErrorVsExternal(
      const CaseResult& cand,
      const ExternalReference& ref,
      const std::vector<double>& output_times,
      double t_skip)
  {
    double worst = 0.0;
    auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
    for (std::size_t i = 0; i < output_times.size() && i < ref.times.size(); ++i)
    {
      if (output_times[i] < t_skip)
        continue;
      if (i >= cand.a_at_output.size())
        break;
      worst = std::max(worst, rel(cand.a_at_output[i], ref.a[i]));
      worst = std::max(worst, rel(cand.c_at_output[i], ref.c[i]));
    }
    return worst;
  }

  // Work-precision sweep at fixed k2: error vs cost across rtol. Accuracy is
  // measured against the external Radau reference (never self-referentially),
  // and the wall-clock repetitions interleave the two methods so thermal or
  // frequency drift cannot bias one side.
  void RunWorkPrecision(double k1, double k2, double k3, const std::vector<double>& output_times, const ExternalReference& ref)
  {
    const std::vector<double> rtols = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
    const int reps = 11;  // odd -> well-defined median
    const double t_skip = 1.0;

    // Cross-solver validation: the tight-tolerance MICM reference used by the
    // stiffness sweep, judged against the independent Radau reference.
    {
      auto tight = RunCase(Method::FullOde, k1, k2, k3, 1e-10, output_times, 0);
      std::cout << "cross-solver agreement (MICM rtol=1e-10 vs Radau rtol=1e-12, all output times): "
                << MaxRelErrorVsExternal(tight, ref, output_times, 0.0) << "\n";
    }

    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    robertson::QssaConstraint constraint(k1, k2, k3);
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

    auto timed_pass = [&](auto& solver, double rtol, bool dae) -> double
    {
      auto state = solver.GetState(1);
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, rtol * 1e-2));
      state.SetCustomRateParameter("r1", k1);
      state.SetCustomRateParameter("r2", k2);
      state.SetCustomRateParameter("r3", k3);
      auto map = state.variable_map_;
      state.variables_[0][map.at("A")] = 1.0;
      state.variables_[0][map.at("B")] = dae ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
      state.variables_[0][map.at("C")] = 0.0;
      state.conditions_[0].temperature_ = 272.5;
      state.conditions_[0].pressure_ = 101253.3;
      state.conditions_[0].air_density_ = 1e6;
      solver.UpdateStateParameters(state);
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
      return std::chrono::duration<double, std::micro>(t1 - t0).count();
    };

    std::ofstream csv("robertson_work_precision.csv");
    csv << "method,rtol,number_of_steps,accepted,rejected,function_calls,"
           "jacobian_updates,decompositions,solves,converged,wallclock_median_us,max_rel_err\n";
    auto write_row = [&](const std::string& method, double rtol, const CaseResult& r, double wall_us, double err)
    {
      csv << method << ',' << rtol << ',' << r.number_of_steps << ',' << r.accepted << ',' << r.rejected << ','
          << r.function_calls << ',' << r.jacobian_updates << ',' << r.decompositions << ',' << r.solves << ','
          << (r.converged ? 1 : 0) << ',' << wall_us << ',' << err << '\n';
    };

    std::cout << "Robertson work-precision at k2=" << k2 << " (errors vs external Radau reference)\n";
    for (double rtol : rtols)
    {
      auto ode = RunCase(Method::FullOde, k1, k2, k3, rtol, output_times, 0);
      auto dae = RunCase(Method::DaeQssa, k1, k2, k3, rtol, output_times, 0);
      std::vector<double> ode_samples, dae_samples;
      for (int rep = 0; rep < reps; ++rep)
      {
        ode_samples.push_back(timed_pass(ode_solver, rtol, false));
        dae_samples.push_back(timed_pass(dae_solver, rtol, true));
      }
      std::sort(ode_samples.begin(), ode_samples.end());
      std::sort(dae_samples.begin(), dae_samples.end());
      const double ode_wall = ode_samples[ode_samples.size() / 2];
      const double dae_wall = dae_samples[dae_samples.size() / 2];
      const double ode_err = MaxRelErrorVsExternal(ode, ref, output_times, t_skip);
      const double dae_err = MaxRelErrorVsExternal(dae, ref, output_times, t_skip);
      write_row("full_ode", rtol, ode, ode_wall, ode_err);
      write_row("dae_qssa", rtol, dae, dae_wall, dae_err);
      std::cout << "rtol=" << rtol << "  ODE err=" << ode_err << " (us=" << ode_wall << ")  DAE err=" << dae_err
                << " (us=" << dae_wall << ")\n";
    }
    std::cout << "wrote robertson_work_precision.csv\n";
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

  const auto reference = LoadExternalReference(std::string(MICM_BENCHMARK_DATA_DIR) + "/robertson_reference.csv");
  if (reference.times.empty())
  {
    std::cout << "no external reference at " << MICM_BENCHMARK_DATA_DIR
              << "/robertson_reference.csv; skipping work-precision sweep "
                 "(generate with benchmark/generate_reference_solutions.py)\n";
  }
  else
  {
    RunWorkPrecision(k1, robertson::K2_DEFAULT, k3, output_times, reference);
  }
  return 0;
}
