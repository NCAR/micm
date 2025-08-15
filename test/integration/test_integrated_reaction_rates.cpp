#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using SparseMatrixTest = micm::SparseMatrix<double>;

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto a = micm::Species("a");
  auto b = micm::Species("b");
  auto c = micm::Species("c");
  auto irr_1 = micm::Species("irr_1");
  auto irr_2 = micm::Species("irr_2");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c, irr_1, irr_2 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1), micm::Yield(irr_1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 =
      micm::ChemicalReactionBuilder()
          .SetReactants({ b })
          .SetProducts({ micm::Yield(c, 1), micm::Yield(irr_2, 1) })
          .SetRateConstant(micm::UserDefinedRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" }))
          .SetPhase(gas_phase)
          .Build();

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                    .SetReactions({ r1, r2 })
                    .Build();

  auto state = solver.GetState(1);

  state.SetCustomRateParameter("r2", 1.0);

  size_t a_idx = state.variable_map_.at(a.name_);
  size_t b_idx = state.variable_map_.at(b.name_);
  size_t c_idx = state.variable_map_.at(c.name_);
  size_t irr1_idx = state.variable_map_.at(irr_1.name_);
  size_t irr2_idx = state.variable_map_.at(irr_2.name_);

  std::vector<double> concentrations(5, 0);

  concentrations[a_idx] = 1.0;

  state.variables_[0] = concentrations;
  state.conditions_[0].temperature_ = 273;
  state.conditions_[0].pressure_ = 1000;

  for (double t{}; t < 100; ++t)
  {
    if (t > 50)
    {
      state.SetCustomRateParameter("r2", 0.0);
    }
    solver.CalculateRateConstants(state);
    std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
    double irr1 = state.variables_[0][irr1_idx];
    double irr2 = state.variables_[0][irr2_idx];
    auto result = solver.Solve(30.0, state);
    EXPECT_GE(state.variables_[0][irr1_idx], irr1);
    EXPECT_GE(state.variables_[0][irr2_idx], irr2);
  }
  std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
}
