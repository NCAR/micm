// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
//
// Test of the TS1 example mechanism (MZ327_TS1.2_20230307 in the Chemistry Cafe)
//
// This only tests that a solver can be built for this mechanism and that it
// can be run for a given set of initial conditions without solver errors.
//
// Comparisons with other chemistry solvers are included in
// https://github.com/NCAR/MUSICA-Performance-Comparison

#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>
#include <micm/solver/rosenbrock.hpp>

TEST(Examples, robertson)
{
  micm::SolverConfig config;
  std::string config_path = "./examples/configs/robertson";
  auto status = config.ReadAndParse(config_path);
  EXPECT_EQ(status, micm::ConfigParseStatus::Success);
  auto solver_params = config.GetSolverParams();
  micm::RosenbrockSolver<> solver{ solver_params.system_,
                                   solver_params.processes_,
                                   micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
  micm::State state = solver.GetState();

  const double temperature = 272.5;  // K
  const double pressure = 101253.3;  // Pa
  const double air_density = 1e6;    // mol m-3

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  std::unordered_map<std::string, std::vector<double>> initial_concentrations = {
    { "A", { 1.0 } },  // mol m-3
    { "B", { 0.0 } },  // mol m-3
    { "C", { 0.0 } },  // mol m-3
  };
  state.SetConcentrations(initial_concentrations);

  std::unordered_map<std::string, std::vector<double>> custom_rate_constants = { { "PHOTO.r1", { 0.04 } },
                                                                                 { "PHOTO.r2", { 3e7 } },
                                                                                 { "PHOTO.r3", { 1e4 } } };

  state.SetCustomRateParameters(custom_rate_constants);

  double time_step = 200.0;  // s

  auto result = solver.Solve(time_step, state);

  EXPECT_EQ(result.state_, (micm::SolverState::Converged));
}
