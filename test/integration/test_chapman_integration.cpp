#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

using SparseMatrixTest = micm::SparseMatrix<double>;

#ifdef USE_YAML
  #include <micm/configure/solver_config.hpp>

TEST(ChapmanIntegration, CanBuildChapmanSystemUsingConfig)
{
  micm::SolverConfig solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  // If parsing fails, it could throw exceptions - we probably want to catch them.
  std::string config_path = "./unit_configs/small_mechanism";
  EXPECT_NO_THROW(solverConfig.ReadAndParse(config_path));

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto options = solver_params.parameters_;  

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(solver_params.system_)
                    .SetReactions(solver_params.processes_)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  micm::State state = solver.GetState();

  EXPECT_EQ(state.relative_tolerance_, 1.0e-6);

  for (size_t n_grid_cell = 0; n_grid_cell < state.number_of_grid_cells_; ++n_grid_cell)
  {
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["Ar"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["CO2"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["H2O"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["N2"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["O1D"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["O"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["O2"]], 1.0e-12);
    EXPECT_EQ(state.absolute_tolerance_[state.variable_map_["O3"]], 1.0e-12);
  }

  // User gives an input of concentrations
  std::unordered_map<std::string, std::vector<double>> concentrations = { { "O", { 0.1 } },   { "O1D", { 0.1 } },
                                                                          { "O2", { 0.1 } },  { "O3", { 0.2 } },
                                                                          { "Ar", { 0.2 } },  { "N2", { 0.3 } },
                                                                          { "H2O", { 0.3 } }, { "CO2", { 0.3 } } };

  state.SetConcentrations(concentrations);

  // User gives an input of photolysis rate constants
  std::unordered_map<std::string, std::vector<double>> photo_rates = { { "PHOTO.O2_1", { 0.1 } },
                                                                       { "PHOTO.O3_1", { 0.2 } },
                                                                       { "PHOTO.O3_2", { 0.3 } } };

  state.SetCustomRateParameters(photo_rates);

  state.conditions_[0].temperature_ = 2;
  state.conditions_[0].pressure_ = 3;
  state.conditions_[0].air_density_ = 0.2;

  for (double t{}; t < 100; ++t)
  {
    state.SetCustomRateParameters(photo_rates);
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(30.0, state);
    // output state
  }
}
#endif

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto o = micm::Species("O");
  auto o1d = micm::Species("O1D");
  auto o2 = micm::Species("O2");
  auto o3 = micm::Species("O3");
  auto m = micm::Species("M");
  auto ar = micm::Species("Ar");
  auto n2 = micm::Species("N2");
  auto h2o = micm::Species("H2O");
  auto co2 = micm::Species("CO2");

  micm::Phase gas_phase{ std::vector<micm::Species>{ o, o1d, o2, o3, m, ar, n2, h2o, co2 } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ o1d, n2 })
                         .SetProducts({ Yields(o, 1), Yields(n2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ o1d, o2 })
                         .SetProducts({ Yields(o, 1), Yields(o2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase);

  micm::Process r3 = micm::Process::Create()
                         .SetReactants({ o, o3 })
                         .SetProducts({ Yields(o2, 2) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 }))
                         .SetPhase(gas_phase);

  micm::Process r4 = micm::Process::Create()
                         .SetReactants({ o, o2, m })
                         .SetProducts({ Yields(o3, 1), Yields(m, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 }))
                         .SetPhase(gas_phase);

  micm::Process photo_1 = micm::Process::Create()
                              .SetReactants({ o2 })
                              .SetProducts({ Yields(o, 2) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO2" }))
                              .SetPhase(gas_phase);

  micm::Process photo_2 = micm::Process::Create()
                              .SetReactants({ o3 })
                              .SetProducts({ Yields(o1d, 1), Yields(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3a" }))
                              .SetPhase(gas_phase);

  micm::Process photo_3 = micm::Process::Create()
                              .SetReactants({ o3 })
                              .SetProducts({ Yields(o, 1), Yields(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3b" }))
                              .SetPhase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ r1, r2, r3, r4, photo_1, photo_2, photo_3 })
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  auto state = solver.GetState();

  std::vector<double> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  state.variables_[0] = concentrations;
  std::vector<double> photo_rates{ 0.1, 0.2, 0.3 };
  state.custom_rate_parameters_[0] = photo_rates;
  state.conditions_[0].temperature_ = 2;
  state.conditions_[0].pressure_ = 3;

  for (double t{}; t < 100; ++t)
  {
    state.custom_rate_parameters_[0] = photo_rates;
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(30.0, state);
    // output state
  }
}
