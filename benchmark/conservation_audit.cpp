// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Element-conservation audit for the tropospheric O3-NOx-HOx benchmark:
// NOy = NO + NO2 + HNO3 is an exact linear invariant of the mechanism (every
// reaction moves nitrogen within that set; N2 is a net-zero spectator), so any
// drift is numerical. QSSA is not conservation-preserving in general — the
// four radical rows are replaced by algebraic balances — but none of the
// constrained radicals carries nitrogen, so the QSSA-DAE should conserve NOy
// as well as the ODE does. This audit quantifies both, across tolerances, over
// the full steady-J integration.
//
// Carbon is deliberately NOT audited: reaction r8 (OH + CH4 -> HO2 + ...)
// drops the methyl carbon by design in this stub mechanism, so carbon closure
// is a property the mechanism itself lacks.
#include "tropospheric_qssa_constraint.hpp"
#include "tropospheric_system.hpp"

#include <micm/CPU.hpp>

#include <algorithm>
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

  enum class Method
  {
    FullOde,
    DaeQssa
  };

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  struct AuditResult
  {
    std::uint64_t accepted = 0;
    bool converged = true;
    double max_noy_drift = 0.0;  // max |NOy(t) - NOy(0)| / NOy(0)
  };

  AuditResult RunAudit(Method method, double rtol, const std::vector<double>& output_times)
  {
    auto sys = tropospheric::MakeSystem();
    const double j1 = tropospheric::JNO2_DEFAULT;
    const double j2 = tropospheric::JO1D_DEFAULT;
    const double j3 = tropospheric::JO3P_DEFAULT;
    auto rates = tropospheric::RadicalRates::At(T_BOX, P_BOX, j1, j2, j3);
    const auto rad0 =
        tropospheric::ProjectRadicals(rates, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O);
    const tropospheric::RefConc ref{ .O1D = rad0.O1D, .O = rad0.O,   .OH = rad0.OH, .HO2 = rad0.HO2,
                                     .O3 = N_O3_0,    .NO = N_NO_0,  .NO2 = N_NO2_0, .CO = N_CO_0,
                                     .CH4 = N_CH4_0,  .N2 = N_N2,    .O2 = N_O2,     .M = N_M,
                                     .H2O = N_H2O };
    tropospheric::QssaRadicalConstraint constraint(rates, ref);

    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = (method == Method::DaeQssa)
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

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
    state.SetCustomRateParameter("p1", j1);
    state.SetCustomRateParameter("p2", j2);
    state.SetCustomRateParameter("p3", j3);
    auto& m = state.variable_map_;
    auto set = [&](const std::string& name, double value) { state.variables_[0][m.at(name)] = value; };
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
    state.conditions_[0].temperature_ = T_BOX;
    state.conditions_[0].pressure_ = P_BOX;
    state.conditions_[0].air_density_ = N_M;
    solver.UpdateStateParameters(state);

    const double noy_0 = N_NO_0 + N_NO2_0;
    AuditResult out;
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
        out.accepted += result.stats_.accepted_;
        done += result.stats_.final_time_;
      }
      current = t_out;
      const double noy =
          state.variables_[0][m.at("NO")] + state.variables_[0][m.at("NO2")] + state.variables_[0][m.at("HNO3")];
      out.max_noy_drift = std::max(out.max_noy_drift, std::abs(noy - noy_0) / noy_0);
    }
    return out;
  }
}  // namespace

int main()
{
  const auto output_times = OutputTimes(1e-3, 1e5, 25);
  const std::vector<double> rtols = { 1e-4, 1e-6, 1e-8 };

  std::ofstream csv("conservation_audit.csv");
  csv << "method,rtol,accepted,converged,max_noy_relative_drift\n";
  csv.precision(12);

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "NOy conservation audit (steady J, j_scale = 1, t -> 1e5 s)\n";
  for (double rtol : rtols)
  {
    auto ode = RunAudit(Method::FullOde, rtol, output_times);
    auto dae = RunAudit(Method::DaeQssa, rtol, output_times);
    csv << "full_ode," << rtol << ',' << ode.accepted << ',' << (ode.converged ? 1 : 0) << ',' << ode.max_noy_drift
        << '\n';
    csv << "dae_qssa," << rtol << ',' << dae.accepted << ',' << (dae.converged ? 1 : 0) << ',' << dae.max_noy_drift
        << '\n';
    std::cout << "rtol=" << rtol << "  ODE NOy drift=" << ode.max_noy_drift << "  DAE NOy drift=" << dae.max_noy_drift
              << (ode.converged && dae.converged ? "" : "  [FAILED]") << "\n";
  }
  std::cout << "wrote conservation_audit.csv\n";
  return 0;
}
