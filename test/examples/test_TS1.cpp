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

TEST(Examples, TS1)
{
  micm::SolverConfig config;
  std::string config_path = "./examples/TS1";
  auto status = config.ReadAndParse(config_path);
  EXPECT_EQ(status, micm::ConfigParseStatus::Success);
  auto solver_params = config.GetSolverParams();
  micm::RosenbrockSolver<> solver{ solver_params.system_,
                                   solver_params.processes_,
                                   micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
  micm::State state = solver.GetState();

  const double temperature = 287.45;                             // K
  const double pressure = 101319.9;                              // Pa
  const double air_density = pressure / (8.3145 * temperature);  // mol m-3

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  // TODO: Replace with CAM-Chem conditions when Derecho is back up
  std::unordered_map<std::string, std::vector<double>> initial_concentrations = {
    { "O2", { 0.21 * air_density } },       // mol m-3
    { "N2", { 0.78 * air_density } },       // mol m-3
    { "HNO3", { 1.0e-9 * air_density } },   // mol m-3
    { "SO2", { 1.0e-9 * air_density } },    // mol m-3
    { "CH4", { 1.8e-6 * air_density } },    // mol m-3
    { "H2O", { 1.0e-3 * air_density } },    // mol m-3
    { "NO2", { 1.0e-9 * air_density } },    // mol m-3
    { "NO", { 1.0e-9 * air_density } },     // mol m-3
    { "OH", { 1.0e-12 * air_density } },    // mol m-3
    { "CL2", { 1.0e-9 * air_density } },    // mol m-3
    { "DMS", { 1.0e-9 * air_density } },    // mol m-3
    { "CO2", { 0.385e-3 * air_density } },  // mol m-3
    { "H2O2", { 1.0e-9 * air_density } },   // mol m-3
    { "ISOP", { 1.0e-9 * air_density } },   // mol m-3
    { "O3", { 50.0e-9 * air_density } },    // mol m-3
    { "CO", { 1.0e-9 * air_density } },     // mol m-3
    { "NH3", { 1.0e-9 * air_density } },    // mol m-3
    { "M", { air_density } }                // mol m-3
  };
  state.SetConcentrations(initial_concentrations);

  // std::cout << "custom rate parameters" << std::endl;
  // for(auto &rc : state.custom_rate_parameter_map_)
  //   std::cout << "    { \"" << rc.first << "\", { 1.0e-8 } }, // s-1" << std::endl;

  // TODO: Replace with CAM-Chem rate constants when Derecho is back up
  std::unordered_map<std::string, std::vector<double>> custom_rate_constants = {
    { "PHOTO.jacet", { 1.0e-8 } },                  // s-1
    { "PHOTO.jalknit->,jch3ooh", { 1.0e-8 } },      // s-1
    { "PHOTO.jalkooh->,jch3ooh", { 1.0e-8 } },      // s-1
    { "PHOTO.jbenzooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jbepomuc->,.10*jno2", { 1.0e-8 } },    // s-1
    { "PHOTO.jbigald->,0.2*jno2", { 1.0e-8 } },     // s-1
    { "PHOTO.jbigald1->,.14*jno2", { 1.0e-8 } },    // s-1
    { "PHOTO.jbigald2->,.20*jno2", { 1.0e-8 } },    // s-1
    { "PHOTO.jbigald3->,.20*jno2", { 1.0e-8 } },    // s-1
    { "PHOTO.jbigald4->,.006*jno2", { 1.0e-8 } },   // s-1
    { "PHOTO.jbrcl", { 1.0e-8 } },                  // s-1
    { "PHOTO.jbro", { 1.0e-8 } },                   // s-1
    { "PHOTO.jbrono2_a", { 1.0e-8 } },              // s-1
    { "PHOTO.jbrono2_b", { 1.0e-8 } },              // s-1
    { "PHOTO.jbzooh->,jch3ooh", { 1.0e-8 } },       // s-1
    { "PHOTO.jc2h5ooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jc3h7ooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jc6h5ooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jccl4", { 1.0e-8 } },                  // s-1
    { "PHOTO.jcf2cl2", { 1.0e-8 } },                // s-1
    { "PHOTO.jcf2clbr", { 1.0e-8 } },               // s-1
    { "PHOTO.jcf3br", { 1.0e-8 } },                 // s-1
    { "PHOTO.jcfc113", { 1.0e-8 } },                // s-1
    { "PHOTO.jcfc114", { 1.0e-8 } },                // s-1
    { "PHOTO.jcfc115", { 1.0e-8 } },                // s-1
    { "PHOTO.jcfcl3", { 1.0e-8 } },                 // s-1
    { "PHOTO.jch2br2", { 1.0e-8 } },                // s-1
    { "PHOTO.jch2o_a", { 1.0e-8 } },                // s-1
    { "PHOTO.jch2o_b", { 1.0e-8 } },                // s-1
    { "PHOTO.jch3br", { 1.0e-8 } },                 // s-1
    { "PHOTO.jch3ccl3", { 1.0e-8 } },               // s-1
    { "PHOTO.jch3cho", { 1.0e-8 } },                // s-1
    { "PHOTO.jch3cl", { 1.0e-8 } },                 // s-1
    { "PHOTO.jch3co3h->,0.28*jh2o2", { 1.0e-8 } },  // s-1
    { "PHOTO.jch3ooh", { 1.0e-8 } },                // s-1
    { "PHOTO.jch4_a", { 1.0e-8 } },                 // s-1
    { "PHOTO.jch4_b", { 1.0e-8 } },                 // s-1
    { "PHOTO.jchbr3", { 1.0e-8 } },                 // s-1
    { "PHOTO.jcl2", { 1.0e-8 } },                   // s-1
    { "PHOTO.jcl2o2", { 1.0e-8 } },                 // s-1
    { "PHOTO.jclo", { 1.0e-8 } },                   // s-1
    { "PHOTO.jclono2_a", { 1.0e-8 } },              // s-1
    { "PHOTO.jclono2_b", { 1.0e-8 } },              // s-1
    { "PHOTO.jco2", { 1.0e-8 } },                   // s-1
    { "PHOTO.jcof2", { 1.0e-8 } },                  // s-1
    { "PHOTO.jcofcl", { 1.0e-8 } },                 // s-1
    { "PHOTO.jeooh->,jch3ooh", { 1.0e-8 } },        // s-1
    { "PHOTO.jglyald", { 1.0e-8 } },                // s-1
    { "PHOTO.jglyoxal->,jmgly", { 1.0e-8 } },       // s-1
    { "PHOTO.jh2402", { 1.0e-8 } },                 // s-1
    { "PHOTO.jh2o2", { 1.0e-8 } },                  // s-1
    { "PHOTO.jh2o_a", { 1.0e-8 } },                 // s-1
    { "PHOTO.jh2o_b", { 1.0e-8 } },                 // s-1
    { "PHOTO.jh2o_c", { 1.0e-8 } },                 // s-1
    { "PHOTO.jh2so4", { 1.0e-8 } },                 // s-1
    { "PHOTO.jhbr", { 1.0e-8 } },                   // s-1
    { "PHOTO.jhcfc141b", { 1.0e-8 } },              // s-1
    { "PHOTO.jhcfc142b", { 1.0e-8 } },              // s-1
    { "PHOTO.jhcfc22", { 1.0e-8 } },                // s-1
    { "PHOTO.jhcl", { 1.0e-8 } },                   // s-1
    { "PHOTO.jhf", { 1.0e-8 } },                    // s-1
    { "PHOTO.jhno3", { 1.0e-8 } },                  // s-1
    { "PHOTO.jho2no2_a", { 1.0e-8 } },              // s-1
    { "PHOTO.jho2no2_b", { 1.0e-8 } },              // s-1
    { "PHOTO.jhobr", { 1.0e-8 } },                  // s-1
    { "PHOTO.jhocl", { 1.0e-8 } },                  // s-1
    { "PHOTO.jhonitr->,jch2o_a", { 1.0e-8 } },      // s-1
    { "PHOTO.jhpald->,.006*jno2", { 1.0e-8 } },     // s-1
    { "PHOTO.jhyac", { 1.0e-8 } },                  // s-1
    { "PHOTO.jisopnooh->,jch3ooh", { 1.0e-8 } },    // s-1
    { "PHOTO.jisopooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jmacr_a", { 1.0e-8 } },                // s-1
    { "PHOTO.jmacr_b", { 1.0e-8 } },                // s-1
    { "PHOTO.jmek->,jacet", { 1.0e-8 } },           // s-1
    { "PHOTO.jmekooh->,jch3ooh", { 1.0e-8 } },      // s-1
    { "PHOTO.jmgly", { 1.0e-8 } },                  // s-1
    { "PHOTO.jmpan->,jpan", { 1.0e-8 } },           // s-1
    { "PHOTO.jmvk", { 1.0e-8 } },                   // s-1
    { "PHOTO.jn2o", { 1.0e-8 } },                   // s-1
    { "PHOTO.jn2o5_a", { 1.0e-8 } },                // s-1
    { "PHOTO.jn2o5_b", { 1.0e-8 } },                // s-1
    { "PHOTO.jnc4cho->,jch2o_a", { 1.0e-8 } },      // s-1
    { "PHOTO.jno2", { 1.0e-8 } },                   // s-1
    { "PHOTO.jno3_a", { 1.0e-8 } },                 // s-1
    { "PHOTO.jno3_b", { 1.0e-8 } },                 // s-1
    { "PHOTO.jno=userdefined,", { 1.0e-8 } },       // s-1
    { "PHOTO.jnoa->,jch2o_a", { 1.0e-8 } },         // s-1
    { "PHOTO.jnterpooh->,jch3ooh", { 1.0e-8 } },    // s-1
    { "PHOTO.jo2_a=userdefined,", { 1.0e-8 } },     // s-1
    { "PHOTO.jo2_b=userdefined,", { 1.0e-8 } },     // s-1
    { "PHOTO.jo3_a", { 1.0e-8 } },                  // s-1
    { "PHOTO.jo3_b", { 1.0e-8 } },                  // s-1
    { "PHOTO.joclo", { 1.0e-8 } },                  // s-1
    { "PHOTO.jocs", { 1.0e-8 } },                   // s-1
    { "PHOTO.jonitr->,jch3cho", { 1.0e-8 } },       // s-1
    { "PHOTO.jpan", { 1.0e-8 } },                   // s-1
    { "PHOTO.jphenooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jpooh->,jch3ooh", { 1.0e-8 } },        // s-1
    { "PHOTO.jrooh->,jch3ooh", { 1.0e-8 } },        // s-1
    { "PHOTO.jsf6", { 1.0e-8 } },                   // s-1
    { "PHOTO.jso", { 1.0e-8 } },                    // s-1
    { "PHOTO.jso2", { 1.0e-8 } },                   // s-1
    { "PHOTO.jso3", { 1.0e-8 } },                   // s-1
    { "PHOTO.jsoa1_a1->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa1_a2->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa2_a1->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa2_a2->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa3_a1->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa3_a2->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa4_a1->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa4_a2->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa5_a1->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jsoa5_a2->,.0004*jno2", { 1.0e-8 } },  // s-1
    { "PHOTO.jtepomuc->,.10*jno2", { 1.0e-8 } },    // s-1
    { "PHOTO.jterp2ooh->,jch3ooh", { 1.0e-8 } },    // s-1
    { "PHOTO.jterpnit->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jterpooh->,jch3ooh", { 1.0e-8 } },     // s-1
    { "PHOTO.jterprd1->,jch3cho", { 1.0e-8 } },     // s-1
    { "PHOTO.jterprd2->,jch3cho", { 1.0e-8 } },     // s-1
    { "PHOTO.jtolooh->,jch3ooh", { 1.0e-8 } },      // s-1
    { "PHOTO.jxooh->,jch3ooh", { 1.0e-8 } },        // s-1
    { "PHOTO.jxylenooh->,jch3ooh", { 1.0e-8 } },    // s-1
    { "PHOTO.jxylolooh->,jch3ooh", { 1.0e-8 } },    // s-1
    { "USER.het1", { 1.0e-8 } },                    // s-1
    { "USER.het2", { 1.0e-8 } },                    // s-1
    { "USER.het3", { 1.0e-8 } },                    // s-1
    { "USER.het4", { 1.0e-8 } },                    // s-1
    { "USER.het5", { 1.0e-8 } },                    // s-1
    { "USER.het6", { 1.0e-8 } },                    // s-1
    { "USER.het7", { 1.0e-8 } },                    // s-1
    { "USER.het8", { 1.0e-8 } },                    // s-1
    { "USER.het9", { 1.0e-8 } },                    // s-1
    { "USER.het10", { 1.0e-8 } },                   // s-1
    { "USER.het11", { 1.0e-8 } },                   // s-1
    { "USER.het12", { 1.0e-8 } },                   // s-1
    { "USER.het13", { 1.0e-8 } },                   // s-1
    { "USER.het14", { 1.0e-8 } },                   // s-1
    { "USER.het15", { 1.0e-8 } },                   // s-1
    { "USER.het16", { 1.0e-8 } },                   // s-1
    { "USER.het17", { 1.0e-8 } },                   // s-1
    { "USER.k_co_oh_jpl19", { 1.0e-8 } }            // s-1
  };
  state.SetCustomRateParameters(custom_rate_constants);

  double time_step = 150.0;  // s

  auto result = solver.Solve(time_step, state);

  EXPECT_EQ(result.state_, (micm::SolverState::Converged));
}