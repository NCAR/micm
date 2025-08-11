#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using SparseMatrixTest = micm::SparseMatrix<double>;

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

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o1d, n2 })
                         .SetProducts({ micm::Yield(o, 1), micm::Yield(n2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o1d, o2 })
                         .SetProducts({ micm::Yield(o, 1), micm::Yield(o2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o, o3 })
                         .SetProducts({ micm::Yield(o2, 2) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o, o2, m })
                         .SetProducts({ micm::Yield(o3, 1), micm::Yield(m, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 }))
                         .SetPhaseName("gas")
                         .Build();

  micm::Process photo_1 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o2 })
                              .SetProducts({ micm::Yield(o, 2) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO2" }))
                              .SetPhaseName("gas")
                              .Build();

  micm::Process photo_2 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o3 })
                              .SetProducts({ micm::Yield(o1d, 1), micm::Yield(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3a" }))
                              .SetPhaseName("gas")
                              .Build();

  micm::Process photo_3 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o3 })
                              .SetProducts({ micm::Yield(o, 1), micm::Yield(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3b" }))
                              .SetPhaseName("gas")
                              .Build();

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ r1, r2, r3, r4, photo_1, photo_2, photo_3 })
                    .SetIgnoreUnusedSpecies(true)
                    .Build();
  ;

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
