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

TEST(Examples, carbon_bond_5)
{
  micm::SolverConfig config;
  std::string config_path = "./examples/configs/carbon_bond_5";
  auto status = config.ReadAndParse(config_path);
  EXPECT_EQ(status, micm::ConfigParseStatus::Success);
  auto solver_params = config.GetSolverParams();
  micm::RosenbrockSolver<> solver{ solver_params.system_,
                                   solver_params.processes_,
                                   micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
  micm::State state = solver.GetState();

  const double temperature = 287.45;  // K
  const double pressure = 101319.9;   // Pa
  const double air_density = 1e6;     // mol m-3

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  std::unordered_map<std::string, std::vector<double>> custom_rate_constants = {
    { "EMIS.NO", { 1.44e-10 } },
    { "EMIS.NO2", { 7.56e-12 } },
    { "EMIS.CO", { 1.9600000000000003e-09 } },
    { "EMIS.SO2", { 1.06e-09 } },
    { "EMIS.FORM", { 1.02e-11 } },
    { "EMIS.MEOH", { 5.920000000000001e-13 } },
    { "EMIS.ALD2", { 4.25e-12 } },
    { "EMIS.PAR", { 4.27e-10 } },
    { "EMIS.ETH", { 4.62e-11 } },
    { "EMIS.OLE", { 1.49e-11 } },
    { "EMIS.IOLE", { 1.49e-11 } },
    { "EMIS.TOL", { 1.53e-11 } },
    { "EMIS.XYL", { 1.4e-11 } },
    { "EMIS.ISOP", { 6.03e-12 } },
    { "PHOTO.NO2", { 0.00477 } },
    { "PHOTO.O3->O1D", { 2.26e-06 } },
    { "PHOTO.O3->O3P", { 0.00025299999999999997 } },
    { "PHOTO.NO3->NO2", { 0.11699999999999999 } },
    { "PHOTO.NO3->NO", { 0.0144 } },
    { "PHOTO.HONO", { 0.000918 } },
    { "PHOTO.H2O2", { 2.59e-06 } },
    { "PHOTO.PNA", { 1.89e-06 } },
    { "PHOTO.HNO3", { 8.61e-08 } },
    { "PHOTO.NTR", { 4.77e-07 } },
    { "PHOTO.ROOH", { 1.81e-06 } },
    { "PHOTO.MEPX", { 1.81e-06 } },
    { "PHOTO.FORM->HO2", { 7.93e-06 } },
    { "PHOTO.FORM->CO", { 2.2e-05 } },
    { "PHOTO.ALD2", { 2.2e-06 } },
    { "PHOTO.PACD", { 1.81e-06 } },
    { "PHOTO.ALDX", { 2.2e-06 } },
    { "PHOTO.OPEN", { 0.0006450000000000001 } },
    { "PHOTO.MGLY", { 7.64e-05 } },
    { "PHOTO.ISPD", { 1.98e-09 } }
  };

  state.SetCustomRateParameters(custom_rate_constants);

  double time_step = 150.0;  // s

  auto result = solver.Solve(time_step, state);

  EXPECT_EQ(result.state_, (micm::SolverState::Converged));
}
