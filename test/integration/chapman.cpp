#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

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

  micm::ArrheniusRateConstantParameters r1_rate_params;
  r1_rate_params.A_ = 2.15e-11;
  r1_rate_params.B_ = 0;
  r1_rate_params.C_ = 110;

  micm::ArrheniusRateConstantParameters r2_rate_params;
  r2_rate_params.A_ = 3.3e-11;
  r2_rate_params.B_ = 0;
  r2_rate_params.C_ = 55;

  micm::ArrheniusRateConstantParameters r3_rate_params;
  r3_rate_params.A_ = 8e-12;
  r3_rate_params.B_ = 0;
  r3_rate_params.C_ = -2060;

  micm::ArrheniusRateConstantParameters r4_rate_params;
  r4_rate_params.A_ = 6.0e-34;
  r4_rate_params.B_ = 0;
  r4_rate_params.C_ = 2.4;

  micm::Process r1 = micm::Process::create()
      .reactants({ o1d, n2 })
      .products({ yields(o, 1), yields(n2, 1) })
      .rate_constant(micm::ArrheniusRateConstant(r1_rate_params))
      .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
      .reactants({ o1d, o2 })
      .products({ yields(o, 1), yields(o2, 1) })
      .rate_constant(micm::ArrheniusRateConstant(r2_rate_params))
      .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
      .reactants({ o, o3 })
      .products({ yields(o2, 2) })
      .rate_constant(micm::ArrheniusRateConstant(r3_rate_params))
      .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
      .reactants({ o, o2, m })
      .products({ yields(o3, 1), yields(m, 1) })
      .rate_constant(micm::ArrheniusRateConstant(r4_rate_params))
      .phase(gas_phase);


  micm::Process photo_1 = micm::Process::create()
                              .reactants({ o2 })
                              .products({ yields(o, 2) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  micm::Process photo_2 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o1d, 1), yields(o2, 1) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  micm::Process photo_3 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o, 1), yields(o2, 1) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  micm::SystemParameters system_params;
  system_params.gas_phase_ = gas_phase;
  micm::RosenbrockSolver solver{ micm::System(system_params),
                                 std::vector<micm::Process>{ r1, r2, r3, r4, photo_1, photo_2, photo_3 },
                                 micm::RosenbrockSolverParameters{} };

  micm::State state = solver.GetState();

  std::vector<double> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  state.variables_[0] = concentrations;
  std::vector<double> photo_rates{ 0.1, 0.2, 0.3 };
  state.custom_rate_parameters_[0] = photo_rates;
  state.conditions_[0].temperature_ = 2;
  state.conditions_[0].pressure_ = 3;

  for (double t{}; t < 100; ++t)
  {
    state.custom_rate_parameters_[0] = photo_rates;
    auto result = solver.Solve(t, t + 0.5, state);
    // output state
  }
}