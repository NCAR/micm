#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

#ifdef USE_JSON
#  include <micm/configure/solver_config.hpp>

TEST(ChapmanIntegration, CanBuildChapmanSystemUsingConfig)
{
  micm::SolverConfig solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  // If parsing fails, it could throw exceptions - we probably want to catch them.
  std::string config_path = "./unit_configs/chapman";
  micm::ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
  EXPECT_EQ(status, micm::ConfigParseStatus::Success);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto options = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters();
  options.ignore_unused_species_ = true;

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{ solver_params.system_,
                                                                 std::move(solver_params.processes_),
                                                                 options };

  micm::State state = solver.GetState();

  // User gives an input of concentrations
  std::unordered_map<std::string, std::vector<double>> concentrations = {
    { "O", { 0.1 } },  { "O1D", { 0.1 } }, { "O2", { 0.1 } },  { "O3", { 0.2 } }, { "M", { 0.2 } },
    { "Ar", { 0.2 } }, { "N2", { 0.3 } },  { "H2O", { 0.3 } }, { "CO2", { 0.3 } }
  };

  state.SetConcentrations(concentrations);

  // User gives an input of photolysis rate constants
  std::unordered_map<std::string, std::vector<double>> photo_rates = { { "PHOTO.O2_1", { 0.1 } },
                                                                       { "PHOTO.O3_1", { 0.2 } },
                                                                       { "PHOTO.O3_2", { 0.3 } } };

  state.SetCustomRateParameters(photo_rates);

  state.conditions_[0].temperature_ = 2;
  state.conditions_[0].pressure_ = 3;

  for (double t{}; t < 100; ++t)
  {
    state.SetCustomRateParameters(photo_rates);
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

  micm::Process r1 = micm::Process::create()
                         .reactants({ o1d, n2 })
                         .products({ yields(o, 1), yields(n2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ o1d, o2 })
                         .products({ yields(o, 1), yields(o2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ o, o3 })
                         .products({ yields(o2, 2) })
                         .rate_constant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ o, o2, m })
                         .products({ yields(o3, 1), yields(m, 1) })
                         .rate_constant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 }))
                         .phase(gas_phase);

  micm::Process photo_1 = micm::Process::create()
                              .reactants({ o2 })
                              .products({ yields(o, 2) })
                              .rate_constant(micm::UserDefinedRateConstant({ .label_ = "jO2" }))
                              .phase(gas_phase);

  micm::Process photo_2 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o1d, 1), yields(o2, 1) })
                              .rate_constant(micm::UserDefinedRateConstant({ .label_ = "jO3a" }))
                              .phase(gas_phase);

  micm::Process photo_3 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o, 1), yields(o2, 1) })
                              .rate_constant(micm::UserDefinedRateConstant({ .label_ = "jO3b" }))
                              .phase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters();
  options.ignore_unused_species_ = true;

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2, r3, r4, photo_1, photo_2, photo_3 },
    options
  };

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
    auto result = solver.Solve(30.0, state);
    // output state
  }
}
