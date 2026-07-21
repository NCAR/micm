// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Diurnal-photolysis benchmark: the tropospheric O3-NOx-HOx mechanism driven
// through a 24-hour day with J(t) updated every 300 s (the outer-loop cadence
// of a host atmosphere model). The QSSA-DAE uses the time-dependent constraint
// (photolysis read from custom rate parameters per evaluation), so a single
// solver instance follows sunrise, noon, sunset, and night.
//
// Three runs share the segmentation: full ODE and QSSA-DAE at rtol 1e-6, and
// a tight full-ODE reference at rtol 1e-10. Per segment the CSV records the
// solar factor, slow-species states, the QSSA-DAE's and working ODE's relative
// deviations from the reference, and the slowest radical timescale
// (tau_HO2 = 1/L_HO2) — the datum a QSSA validity monitor would compare
// against the step size. The summary reports the DAE error envelope by
// day / twilight / night regime.
#include "diurnal_qssa_constraint.hpp"
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

  constexpr double DAY_SECONDS = 86400.0;
  constexpr double SEGMENT_SECONDS = 300.0;
  constexpr double SUNRISE = 6.0 * 3600.0;
  constexpr double SUNSET = 18.0 * 3600.0;

  /// Solar scale factor: half-sine daylight between 06:00 and 18:00, zero at night.
  double SolarFactor(double t)
  {
    if (t < SUNRISE || t > SUNSET)
      return 0.0;
    return std::sin(M_PI * (t - SUNRISE) / (SUNSET - SUNRISE));
  }

  enum class Method
  {
    FullOde,
    DaeQssa
  };

  struct Trajectory
  {
    std::vector<double> o3, no, no2, hno3, oh, ho2;
    bool converged = true;
  };

  // `project_radicals` is the DAE's operational recipe: before each segment,
  // radicals are warm-started onto the QSSA manifold analytically (damped
  // fixed point). This is not optional at sunrise — overnight NO titration
  // (NO + O3 -> NO2) leaves NO ~ 0, so at HO2 = 0 the HO2 constraint row's
  // Jacobian -(k_HO2_NO*NO + 4*k_HO2_HO2*HO2) vanishes and the initialization
  // Newton sits exactly on the documented z = 0 square-root singularity; the
  // analytic projection handles that manifold branch in closed form.
  template<typename Solver>
  Trajectory RunDay(Solver& solver, double rtol, int segments, bool project_radicals)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
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
    // Midnight start: J = 0, radicals on the (trivial) night manifold.
    set("O1D", 0.0);
    set("O", 0.0);
    set("OH", 0.0);
    set("HO2", 0.0);
    state.conditions_[0].temperature_ = T_BOX;
    state.conditions_[0].pressure_ = P_BOX;
    state.conditions_[0].air_density_ = N_M;
    solver.UpdateStateParameters(state);

    Trajectory out;
    for (int segment = 0; segment < segments; ++segment)
    {
      const double t_mid = (segment + 0.5) * SEGMENT_SECONDS;
      const double s = SolarFactor(t_mid);
      state.SetCustomRateParameter("p1", s * tropospheric::JNO2_DEFAULT);
      state.SetCustomRateParameter("p2", s * tropospheric::JO1D_DEFAULT);
      state.SetCustomRateParameter("p3", s * tropospheric::JO3P_DEFAULT);
      if (project_radicals)
      {
        if (s <= 0.0)
        {
          set("O1D", 0.0);
          set("O", 0.0);
          set("OH", 0.0);
          set("HO2", 0.0);
        }
        else
        {
          auto rates_now = tropospheric::RadicalRates::At(
              T_BOX, P_BOX, s * tropospheric::JNO2_DEFAULT, s * tropospheric::JO1D_DEFAULT, s * tropospheric::JO3P_DEFAULT);
          auto c = tropospheric::ProjectRadicals(
              rates_now,
              state.variables_[0][m.at("O3")],
              state.variables_[0][m.at("NO")],
              state.variables_[0][m.at("NO2")],
              state.variables_[0][m.at("CO")],
              state.variables_[0][m.at("CH4")],
              state.variables_[0][m.at("N2")],
              state.variables_[0][m.at("O2")],
              state.variables_[0][m.at("M")],
              state.variables_[0][m.at("H2O")]);
          set("O1D", c.O1D);
          set("O", c.O);
          set("OH", c.OH);
          set("HO2", c.HO2);
        }
      }
      solver.UpdateStateParameters(state);
      double done = 0.0;
      int guard = 0;
      while (done < SEGMENT_SECONDS && guard++ < 100000)
      {
        auto result = solver.Solve(SEGMENT_SECONDS - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
        {
          std::cout << "segment " << segment << " (t=" << segment * SEGMENT_SECONDS << " s, s=" << s
                    << "): " << micm::SolverStateToString(result.state_) << "\n";
          out.converged = false;
          return out;
        }
        if (result.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          return out;
        }
        done += result.stats_.final_time_;
      }
      out.o3.push_back(state.variables_[0][m.at("O3")]);
      out.no.push_back(state.variables_[0][m.at("NO")]);
      out.no2.push_back(state.variables_[0][m.at("NO2")]);
      out.hno3.push_back(state.variables_[0][m.at("HNO3")]);
      out.oh.push_back(state.variables_[0][m.at("OH")]);
      out.ho2.push_back(state.variables_[0][m.at("HO2")]);
    }
    return out;
  }

  /// The guarded DAE recipe: QSSA-DAE while the sun is up, full-ODE fallback
  /// at night. Once overnight titration removes NO, the HO2 balance collapses
  /// to 0 = -2 k HO2^2 — a degenerate double root whose Jacobian is singular
  /// at the solution — so night QSSA is structurally invalid (the validity
  /// monitor's tau_HO2 -> infinity signals exactly this). Variables are copied
  /// between the two solvers' states at the sunrise/sunset switches.
  template<typename OdeSolver, typename DaeSolver>
  Trajectory RunGuardedDay(OdeSolver& ode_solver, DaeSolver& dae_solver, double rtol, int segments)
  {
    auto ode_state = ode_solver.GetState(1);
    auto dae_state = dae_solver.GetState(1);
    for (auto* state_ptr : { static_cast<void*>(&ode_state), static_cast<void*>(&dae_state) })
    {
      (void)state_ptr;
    }
    auto setup = [&](auto& state)
    {
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
      auto& mm = state.variable_map_;
      auto set = [&](const std::string& name, double value) { state.variables_[0][mm.at(name)] = value; };
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
      set("O1D", 0.0);
      set("O", 0.0);
      set("OH", 0.0);
      set("HO2", 0.0);
      state.conditions_[0].temperature_ = T_BOX;
      state.conditions_[0].pressure_ = P_BOX;
      state.conditions_[0].air_density_ = N_M;
    };
    setup(ode_state);
    setup(dae_state);

    Trajectory out;
    bool using_dae = false;
    for (int segment = 0; segment < segments; ++segment)
    {
      const double t_mid = (segment + 0.5) * SEGMENT_SECONDS;
      const double s = SolarFactor(t_mid);
      const bool want_dae = s > 0.0;
      if (want_dae != using_dae)
      {
        // Hand the 16 species across at the regime switch.
        auto& from = using_dae ? dae_state.variables_ : ode_state.variables_;
        auto& from_map = using_dae ? dae_state.variable_map_ : ode_state.variable_map_;
        auto& to = want_dae ? dae_state.variables_ : ode_state.variables_;
        auto& to_map = want_dae ? dae_state.variable_map_ : ode_state.variable_map_;
        for (const auto& [name, from_index] : from_map)
        {
          to[0][to_map.at(name)] = from[0][from_index];
        }
        using_dae = want_dae;
      }

      auto run_segment = [&](auto& solver, auto& state) -> bool
      {
        auto& mm = state.variable_map_;
        auto set = [&](const std::string& name, double value) { state.variables_[0][mm.at(name)] = value; };
        state.SetCustomRateParameter("p1", s * tropospheric::JNO2_DEFAULT);
        state.SetCustomRateParameter("p2", s * tropospheric::JO1D_DEFAULT);
        state.SetCustomRateParameter("p3", s * tropospheric::JO3P_DEFAULT);
        if (want_dae)
        {
          auto rates_now = tropospheric::RadicalRates::At(
              T_BOX, P_BOX, s * tropospheric::JNO2_DEFAULT, s * tropospheric::JO1D_DEFAULT, s * tropospheric::JO3P_DEFAULT);
          auto c = tropospheric::ProjectRadicals(
              rates_now,
              state.variables_[0][mm.at("O3")],
              state.variables_[0][mm.at("NO")],
              state.variables_[0][mm.at("NO2")],
              state.variables_[0][mm.at("CO")],
              state.variables_[0][mm.at("CH4")],
              state.variables_[0][mm.at("N2")],
              state.variables_[0][mm.at("O2")],
              state.variables_[0][mm.at("M")],
              state.variables_[0][mm.at("H2O")]);
          set("O1D", c.O1D);
          set("O", c.O);
          set("OH", c.OH);
          set("HO2", c.HO2);
        }
        solver.UpdateStateParameters(state);
        double done = 0.0;
        int guard = 0;
        while (done < SEGMENT_SECONDS && guard++ < 100000)
        {
          auto result = solver.Solve(SEGMENT_SECONDS - done, state);
          if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
          {
            std::cout << "guarded segment " << segment << " (s=" << s
                      << "): " << micm::SolverStateToString(result.state_) << "\n";
            return false;
          }
          if (result.stats_.final_time_ <= 0.0)
            return false;
          done += result.stats_.final_time_;
        }
        auto record = [&](const std::string& name) { return state.variables_[0][mm.at(name)]; };
        out.o3.push_back(record("O3"));
        out.no.push_back(record("NO"));
        out.no2.push_back(record("NO2"));
        out.hno3.push_back(record("HNO3"));
        out.oh.push_back(record("OH"));
        out.ho2.push_back(record("HO2"));
        return true;
      };

      const bool ok = want_dae ? run_segment(dae_solver, dae_state) : run_segment(ode_solver, ode_state);
      if (!ok)
      {
        out.converged = false;
        return out;
      }
    }
    return out;
  }
}  // namespace

int main()
{
  const int segments = static_cast<int>(DAY_SECONDS / SEGMENT_SECONDS);
  const double rtol = 1e-6;

  auto sys = tropospheric::MakeSystem();
  // Reference rates and per-row scaling from the noon radical state.
  auto noon_rates = tropospheric::RadicalRates::At(
      T_BOX, P_BOX, tropospheric::JNO2_DEFAULT, tropospheric::JO1D_DEFAULT, tropospheric::JO3P_DEFAULT);
  const auto noon_rad =
      tropospheric::ProjectRadicals(noon_rates, N_O3_0, N_NO_0, N_NO2_0, N_CO_0, N_CH4_0, N_N2, N_O2, N_M, N_H2O);
  const tropospheric::RefConc ref{ .O1D = noon_rad.O1D, .O = noon_rad.O,   .OH = noon_rad.OH, .HO2 = noon_rad.HO2,
                                   .O3 = N_O3_0,        .NO = N_NO_0,      .NO2 = N_NO2_0,    .CO = N_CO_0,
                                   .CH4 = N_CH4_0,      .N2 = N_N2,        .O2 = N_O2,        .M = N_M,
                                   .H2O = N_H2O };
  tropospheric::DiurnalQssaRadicalConstraint constraint(noon_rates, ref);

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

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "Diurnal tropospheric benchmark: 24 h, J updated every " << SEGMENT_SECONDS << " s\n";
  auto reference = RunDay(ode_solver, 1e-10, segments, false);
  auto ode = RunDay(ode_solver, rtol, segments, false);
  auto dae = RunGuardedDay(ode_solver, dae_solver, rtol, segments);
  if (!reference.converged || !ode.converged || !dae.converged)
  {
    std::cout << "FAILED: reference=" << reference.converged << " ode=" << ode.converged << " dae=" << dae.converged
              << "\n";
    return 1;
  }

  std::ofstream csv("diurnal_dae.csv");
  csv << std::scientific << std::setprecision(12);
  csv << "hour,solar_factor,o3_ref,no_ref,no2_ref,hno3_ref,oh_dae,ho2_dae,tau_ho2_s,"
         "ode_max_rel_err,dae_max_rel_err\n";

  auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
  double day_worst = 0.0, twilight_worst = 0.0, night_worst = 0.0;
  double ode_worst = 0.0;
  for (int segment = 0; segment < segments; ++segment)
  {
    const double t_end = (segment + 1) * SEGMENT_SECONDS;
    const double s = SolarFactor(t_end - 0.5 * SEGMENT_SECONDS);
    // Error over the slow observables (O3, NO2, HNO3 once formed).
    auto seg_err = [&](const Trajectory& traj)
    {
      double worst = rel(traj.o3[segment], reference.o3[segment]);
      worst = std::max(worst, rel(traj.no2[segment], reference.no2[segment]));
      if (reference.hno3[segment] > 1e6)
        worst = std::max(worst, rel(traj.hno3[segment], reference.hno3[segment]));
      return worst;
    };
    const double dae_err = seg_err(dae);
    const double ode_err = seg_err(ode);
    ode_worst = std::max(ode_worst, ode_err);
    if (s <= 0.0)
      night_worst = std::max(night_worst, dae_err);
    else if (s < 0.2)
      twilight_worst = std::max(twilight_worst, dae_err);
    else
      day_worst = std::max(day_worst, dae_err);

    const double L_HO2 = noon_rates.k_HO2_NO * dae.no[segment] + 2.0 * noon_rates.k_HO2_HO2 * dae.ho2[segment];
    const double tau_ho2 = L_HO2 > 0.0 ? 1.0 / L_HO2 : std::numeric_limits<double>::infinity();
    csv << t_end / 3600.0 << ',' << s << ',' << reference.o3[segment] << ',' << reference.no[segment] << ','
        << reference.no2[segment] << ',' << reference.hno3[segment] << ',' << dae.oh[segment] << ',' << dae.ho2[segment]
        << ',' << tau_ho2 << ',' << ode_err << ',' << dae_err << '\n';
  }

  std::cout << "QSSA-DAE error envelope vs tight-ODE reference (slow observables):\n";
  std::cout << "  night (s=0):        " << night_worst << "\n";
  std::cout << "  twilight (0<s<0.2): " << twilight_worst << "\n";
  std::cout << "  day (s>=0.2):       " << day_worst << "\n";
  std::cout << "  working-ODE envelope (all regimes): " << ode_worst << "\n";
  std::cout << "wrote diurnal_dae.csv\n";
  return 0;
}
