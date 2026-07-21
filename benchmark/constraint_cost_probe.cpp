// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Phase-3 profiling probe: where does the QSSA-DAE's per-step wall-clock
// overhead come from? Three variants of the Robertson integration are timed
// with interleaved repetitions:
//
//   ode        - plain full ODE (baseline)
//   ode+null   - full ODE plus a null external model whose ODE-side forcing
//                and Jacobian callbacks do nothing: isolates the pure
//                std::function dispatch overhead added to the ProcessSet
//                sweeps, with bitwise-identical dynamics (asserted via equal
//                accepted-step counts)
//   dae        - the QSSA-DAE (constraint math + mass coupling + per-call
//                constraint initialization on top of dispatch)
//
// Deterministic counters are printed alongside so the extra constraint-
// initialization factorizations (decompositions - steps) are visible.
#include "robertson_qssa_constraint.hpp"
#include "robertson_system.hpp"

#include <micm/CPU.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  /// ODE-side external model whose callbacks do nothing: measures the cost of
  /// routing the ProcessSet sweeps through external-model std::functions.
  class NullExternalModel
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "A" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>&,
        const SparseMatrixPolicy&) const
    {
      return [](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&) { };
    }
  };

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  struct Counters
  {
    std::uint64_t accepted = 0;
    std::uint64_t function_calls = 0;
    std::uint64_t jacobian_updates = 0;
    std::uint64_t decompositions = 0;
    std::uint64_t solves = 0;
    std::uint64_t constraint_init_iterations = 0;
  };

  template<typename Solver>
  void InitState(Solver& solver, auto& state, double k1, double k2, double k3, double rtol, bool dae)
  {
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, rtol * 1e-2));
    state.SetCustomRateParameter("r1", k1);
    state.SetCustomRateParameter("r2", k2);
    state.SetCustomRateParameter("r3", k3);
    const auto& map = state.variable_map_;
    state.variables_[0][map.at("A")] = 1.0;
    state.variables_[0][map.at("B")] = dae ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
    state.variables_[0][map.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 272.5;
    state.conditions_[0].pressure_ = 101253.3;
    state.conditions_[0].air_density_ = 1e6;
    solver.UpdateStateParameters(state);
  }

  template<typename Solver>
  Counters CountersPass(Solver& solver, double k1, double k2, double k3, double rtol, bool dae, const std::vector<double>& output_times)
  {
    auto state = solver.GetState(1);
    InitState(solver, state, k1, k2, k3, rtol, dae);
    Counters out;
    double current = 0.0;
    for (double t_out : output_times)
    {
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
          return out;
        if (result.stats_.final_time_ <= 0.0)
          return out;
        out.accepted += result.stats_.accepted_;
        out.function_calls += result.stats_.function_calls_;
        out.jacobian_updates += result.stats_.jacobian_updates_;
        out.decompositions += result.stats_.decompositions_;
        out.solves += result.stats_.solves_;
        out.constraint_init_iterations += result.stats_.constraint_init_iterations_;
        done += result.stats_.final_time_;
      }
      current = t_out;
    }
    return out;
  }

  template<typename Solver>
  double TimedPass(Solver& solver, double k1, double k2, double k3, double rtol, bool dae, const std::vector<double>& output_times)
  {
    auto state = solver.GetState(1);
    InitState(solver, state, k1, k2, k3, rtol, dae);
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
  }
}  // namespace

int main()
{
  const double k1 = robertson::K1_DEFAULT;
  const double k2 = robertson::K2_DEFAULT;
  const double k3 = robertson::K3_DEFAULT;
  const double rtol = 1e-6;
  const int reps = 21;
  const auto output_times = OutputTimes(1e-3, 1e6, 19);

  auto sys = robertson::MakeSystem();
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  robertson::QssaConstraint constraint(k1, k2, k3);

  auto ode = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                 .SetSystem(micm::System(sys.gas_phase))
                 .SetReactions(sys.processes)
                 .SetReorderState(false)
                 .Build();
  auto ode_null = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                      .SetSystem(micm::System(sys.gas_phase))
                      .SetReactions(sys.processes)
                      .AddExternalModel(NullExternalModel())
                      .SetReorderState(false)
                      .Build();
  auto dae = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                 .SetSystem(micm::System(sys.gas_phase))
                 .SetReactions(sys.processes)
                 .AddExternalModel(constraint)
                 .SetReorderState(false)
                 .Build();

  const auto ode_counters = CountersPass(ode, k1, k2, k3, rtol, false, output_times);
  const auto null_counters = CountersPass(ode_null, k1, k2, k3, rtol, false, output_times);
  const auto dae_counters = CountersPass(dae, k1, k2, k3, rtol, true, output_times);

  if (ode_counters.accepted != null_counters.accepted)
  {
    std::cerr << "null external model changed the integration (accepted " << ode_counters.accepted << " vs "
              << null_counters.accepted << ") — dispatch measurement invalid\n";
    return 1;
  }

  std::vector<double> t_ode, t_null, t_dae;
  for (int rep = 0; rep < reps; ++rep)
  {
    t_ode.push_back(TimedPass(ode, k1, k2, k3, rtol, false, output_times));
    t_null.push_back(TimedPass(ode_null, k1, k2, k3, rtol, false, output_times));
    t_dae.push_back(TimedPass(dae, k1, k2, k3, rtol, true, output_times));
  }
  auto median = [](std::vector<double>& v)
  {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  const double us_ode = median(t_ode), us_null = median(t_null), us_dae = median(t_dae);

  auto report = [&](const char* name, const Counters& c, double us)
  {
    std::cout << name << ": us=" << us << "  accepted=" << c.accepted << "  f_calls=" << c.function_calls
              << "  jac=" << c.jacobian_updates << "  decomp=" << c.decompositions << "  solves=" << c.solves
              << "  init_iters=" << c.constraint_init_iterations << "  us/step=" << us / c.accepted << "\n";
  };
  std::cout << "Robertson constraint-cost probe (rtol=" << rtol << ", k2=" << k2 << ", median of " << reps
            << " interleaved reps)\n";
  report("ode     ", ode_counters, us_ode);
  report("ode+null", null_counters, us_null);
  report("dae     ", dae_counters, us_dae);

  const double dispatch_us = us_null - us_ode;
  const double dae_extra_us = us_dae - us_ode;
  std::cout << "\ndecomposition of the DAE overhead (" << dae_extra_us << " us total, "
            << 100.0 * dae_extra_us / us_ode << "% of baseline):\n";
  std::cout << "  external-model dispatch (null model):        " << dispatch_us << " us\n";
  std::cout << "  constraint math + mass coupling + init:      " << dae_extra_us - dispatch_us << " us\n";
  std::cout << "  extra factorizations vs steps (init Newton): "
            << (dae_counters.decompositions - dae_counters.accepted) << " of " << dae_counters.decompositions
            << " total decompositions\n";

  // Single-call variant: one Solve(1e6) instead of 19 dense-output segments,
  // so constraint initialization runs once and the per-step machinery cost
  // (constraint residual/Jacobian sweeps + mass coupling) stands alone.
  const std::vector<double> single_time = { 1e6 };
  const auto ode_single = CountersPass(ode, k1, k2, k3, rtol, false, single_time);
  const auto dae_single = CountersPass(dae, k1, k2, k3, rtol, true, single_time);
  std::vector<double> t_ode_single, t_dae_single;
  for (int rep = 0; rep < reps; ++rep)
  {
    t_ode_single.push_back(TimedPass(ode, k1, k2, k3, rtol, false, single_time));
    t_dae_single.push_back(TimedPass(dae, k1, k2, k3, rtol, true, single_time));
  }
  const double us_ode_single = median(t_ode_single);
  const double us_dae_single = median(t_dae_single);
  const double ode_per_step = us_ode_single / ode_single.accepted;
  const double dae_per_step = us_dae_single / dae_single.accepted;
  std::cout << "\nsingle-call variant (init amortized to one call):\n";
  report("ode     ", ode_single, us_ode_single);
  report("dae     ", dae_single, us_dae_single);
  std::cout << "  per-step DAE machinery overhead: " << dae_per_step - ode_per_step << " us/step ("
            << 100.0 * (dae_per_step - ode_per_step) / ode_per_step << "% of baseline per-step)\n";
  std::cout << "  implied per-init cost in the dense-output run: "
            << (dae_extra_us - dispatch_us - (dae_per_step - ode_per_step) * dae_counters.accepted) / 19.0
            << " us per Solve() call\n";
  return 0;
}
