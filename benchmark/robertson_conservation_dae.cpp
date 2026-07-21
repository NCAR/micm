// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Robertson in its conservation-DAE form (IVP test-set style, t -> 1e11):
// the C row is replaced by the mass invariant A + B + C = 1 using MICM's
// built-in LinearConstraint, with C algebraic and A, B differential.
//
// Unlike the QSSA benchmark, this DAE is EXACT — conservation is an invariant
// of the ODE — so both formulations converge to the same solution and the
// comparison isolates solver behavior over fourteen decades of time:
// terminal accuracy against the external Radau reference (which reproduces
// the published IVP test-set values to ~8 digits) and drift of the mass
// invariant, which the full ODE only conserves to accumulated truncation
// error while the DAE enforces it by construction.
#include "robertson_system.hpp"

#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/linear_constraint.hpp>

#include <algorithm>
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
    ConservationDae
  };

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
    std::vector<double> a_at_output;
    std::vector<double> b_at_output;
    std::vector<double> c_at_output;
  };

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  CaseResult RunCase(Method method, double k1, double k2, double k3, double rtol, const std::vector<double>& output_times)
  {
    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

    // The constraint names the species objects, so they must come from the
    // same system build the solver uses.
    micm::Species a_species("A"), b_species("B"), c_species("C");
    std::vector<micm::Constraint> constraints;
    if (method == Method::ConservationDae)
    {
      constraints.push_back(micm::LinearConstraint(
          "mass", c_species, { { a_species, 1.0 }, { b_species, 1.0 }, { c_species, 1.0 } }, 1.0));
    }

    auto solver = (method == Method::ConservationDae)
                      ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(micm::System(sys.gas_phase))
                            .SetReactions(sys.processes)
                            .SetConstraints(std::move(constraints))
                            .SetReorderState(false)
                            .Build()
                      : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                            .SetSystem(micm::System(sys.gas_phase))
                            .SetReactions(sys.processes)
                            .SetReorderState(false)
                            .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    // Per-species absolute floors: B falls to ~8e-14 by t = 1e11 and needs a
    // floor far below that to stay relatively resolved.
    {
      const auto& map = state.variable_map_;
      std::vector<double> atol(state.state_size_, 1.0e-14);
      atol[map.at("B")] = 1.0e-18;
      state.SetAbsoluteTolerances(atol);
    }
    state.SetCustomRateParameter("r1", k1);
    state.SetCustomRateParameter("r2", k2);
    state.SetCustomRateParameter("r3", k3);
    const auto& map = state.variable_map_;
    state.variables_[0][map.at("A")] = 1.0;
    state.variables_[0][map.at("B")] = 0.0;
    state.variables_[0][map.at("C")] = 0.0;  // consistent: A + B + C = 1
    state.conditions_[0].temperature_ = 272.5;
    state.conditions_[0].pressure_ = 101253.3;
    state.conditions_[0].air_density_ = 1e6;
    solver.UpdateStateParameters(state);

    CaseResult out;
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
      out.a_at_output.push_back(state.variables_[0][map.at("A")]);
      out.b_at_output.push_back(state.variables_[0][map.at("B")]);
      out.c_at_output.push_back(state.variables_[0][map.at("C")]);
      current = t_out;
      if (!out.converged)
        break;
    }
    return out;
  }

  struct ExternalReference
  {
    std::vector<double> times;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
  };

  ExternalReference LoadExternalReference(const std::string& path)
  {
    ExternalReference ref;
    std::ifstream file(path);
    if (!file)
      return ref;
    std::string line;
    std::getline(file, line);  // header: time,A,B,C
    while (std::getline(file, line))
    {
      std::istringstream row(line);
      std::string field;
      std::vector<double> values;
      while (std::getline(row, field, ','))
        values.push_back(std::stod(field));
      if (values.size() < 4)
        continue;
      ref.times.push_back(values[0]);
      ref.a.push_back(values[1]);
      ref.b.push_back(values[2]);
      ref.c.push_back(values[3]);
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
      worst = std::max(worst, rel(cand.b_at_output[i], ref.b[i]));
      worst = std::max(worst, rel(cand.c_at_output[i], ref.c[i]));
    }
    return worst;
  }

  double MaxConservationDrift(const CaseResult& cand)
  {
    double worst = 0.0;
    for (std::size_t i = 0; i < cand.a_at_output.size(); ++i)
      worst = std::max(worst, std::abs(cand.a_at_output[i] + cand.b_at_output[i] + cand.c_at_output[i] - 1.0));
    return worst;
  }
}  // namespace

int main()
{
  const double k1 = robertson::K1_DEFAULT;
  const double k2 = robertson::K2_DEFAULT;
  const double k3 = robertson::K3_DEFAULT;
  const double t_skip = 1.0;
  const auto output_times = OutputTimes(1e-3, 1e11, 29);
  const std::vector<double> rtols = { 1e-4, 1e-6, 1e-8, 1e-10 };

  const auto reference =
      LoadExternalReference(std::string(MICM_BENCHMARK_DATA_DIR) + "/robertson_longtime_reference.csv");
  if (reference.times.empty())
  {
    std::cout << "no external reference at " << MICM_BENCHMARK_DATA_DIR
              << "/robertson_longtime_reference.csv; generate with benchmark/generate_reference_solutions.py\n";
    return 1;
  }

  std::ofstream csv("robertson_conservation.csv");
  csv << "method,rtol,number_of_steps,accepted,rejected,function_calls,"
         "jacobian_updates,decompositions,solves,converged,terminal_a,terminal_b,terminal_c,"
         "max_conservation_drift,max_rel_err\n";
  csv.precision(12);

  auto write_row = [&](const std::string& method, double rtol, const CaseResult& r, double drift, double err)
  {
    csv << method << ',' << rtol << ',' << r.number_of_steps << ',' << r.accepted << ',' << r.rejected << ','
        << r.function_calls << ',' << r.jacobian_updates << ',' << r.decompositions << ',' << r.solves << ','
        << (r.converged ? 1 : 0) << ',' << r.a_at_output.back() << ',' << r.b_at_output.back() << ','
        << r.c_at_output.back() << ',' << drift << ',' << err << '\n';
  };

  std::cout << "Robertson conservation-DAE (t -> 1e11, errors vs external Radau reference,\n"
               "which matches the published IVP test-set terminal values to ~8 digits)\n";
  for (double rtol : rtols)
  {
    auto ode = RunCase(Method::FullOde, k1, k2, k3, rtol, output_times);
    auto dae = RunCase(Method::ConservationDae, k1, k2, k3, rtol, output_times);
    const double ode_drift = MaxConservationDrift(ode);
    const double dae_drift = MaxConservationDrift(dae);
    const double ode_err = MaxRelErrorVsExternal(ode, reference, output_times, t_skip);
    const double dae_err = MaxRelErrorVsExternal(dae, reference, output_times, t_skip);
    write_row("full_ode", rtol, ode, ode_drift, ode_err);
    write_row("conservation_dae", rtol, dae, dae_drift, dae_err);
    std::cout << "rtol=" << rtol << "  ODE steps=" << ode.accepted << " err=" << ode_err << " drift=" << ode_drift
              << "  DAE steps=" << dae.accepted << " err=" << dae_err << " drift=" << dae_drift
              << (ode.converged && dae.converged ? "" : "  [FAILED]") << "\n";
  }
  std::cout << "wrote robertson_conservation.csv\n";
  return 0;
}
