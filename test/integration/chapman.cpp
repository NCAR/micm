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

  micm::Process r1{ .reactants_ = { o1d, n2 },
                    .products_ = { yields(o, 1), yields(n2, 1) },
                    .rate_constant_ =
                        micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .C_ = 110 }),
                    .phase_ = &gas_phase };

  micm::Process r2{ .reactants_ = { o1d, o2 },
                    .products_ = { yields(o, 1), yields(o2, 1) },
                    .rate_constant_ =
                        micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .C_ = 55 }),
                    .phase_ = &gas_phase };

  micm::Process r3{ .reactants_ = { o, o3 },
                    .products_ = { yields(o2, 2) },
                    .rate_constant_ =
                        micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .C_ = -2060 }),
                    .phase_ = &gas_phase };

  micm::Process r4{ .reactants_ = { o, o2, m },
                    .products_ = { yields(o3, 1), yields(m, 1) },
                    .rate_constant_ =
                        micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .C_ = 2.4 }),
                    .phase_ = &gas_phase };

  micm::Process photo_1{ .reactants_ = { o2 },
                         .products_ = { yields(o, 2) },
                         .rate_constant_ = micm::PhotolysisRateConstant(),
                         .phase_ = &gas_phase };

  micm::Process photo_2{ .reactants_ = { o3 },
                         .products_ = { yields(o1d, 1), yields(o2, 1) },
                         .rate_constant_ = micm::PhotolysisRateConstant(),
                         .phase_ = &gas_phase };

  micm::Process photo_3{ .reactants_ = { o3 },
                         .products_ = { yields(o, 1), yields(o2, 1) },
                         .rate_constant_ = micm::PhotolysisRateConstant(),
                         .phase_ = &gas_phase };

  micm::State state{
    micm::System(gas_phase),

    std::vector<micm::Process>{ r1, r2, r3, r4, photo_1, photo_2, photo_3 },
  };

  state.set_concentrations({ { "O", 0.1 },
                             { "O1D", 0.1 },
                             { "O2", 0.1 },
                             { "O3", 0.2 },
                             { "M", 0.2 },
                             { "Ar", 0.2 },
                             { "N2", 0.3 },
                             { "H2O", 0.3 },
                             { "CO2", 0.3 } });
  state.temperature_ = 2;
  state.pressure_ = 3;

  double photo_rates[3] = { 0.1, 0.2, 0.3 };

  micm::RosenbrockSolver solver(state);
  state.update_photo_rates(photo_rates);

  for (double t{}; t < 100; ++t)
  {
    state.update_photo_rates(photo_rates);
    auto result = solver.Solve(t, t + 0.5, state);
    // output state
  }
}