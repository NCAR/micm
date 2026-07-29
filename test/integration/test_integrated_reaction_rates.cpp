#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using SparseMatrixTest = micm::SparseMatrix<micm::Real>;

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto a = micm::Species("a");
  auto b = micm::Species("b");
  auto c = micm::Species("c");
  auto irr_1 = micm::Species("irr_1");
  auto irr_2 = micm::Species("irr_2");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c, irr_1, irr_2 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::StoichSpecies(b, 1), micm::StoichSpecies(irr_1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::StoichSpecies(c, 1), micm::StoichSpecies(irr_2, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" })
                         .SetPhase(gas_phase)
                         .Build();

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(micm::System(gas_phase))
                    .SetReactions({ r1, r2 })
                    .Build();

  auto state = solver.GetState(1);

  state.SetCustomRateParameter("r2", 1.0);

  micm::Index a_idx = state.variable_map_.at(a.name_);
  micm::Index b_idx = state.variable_map_.at(b.name_);
  micm::Index c_idx = state.variable_map_.at(c.name_);
  micm::Index irr1_idx = state.variable_map_.at(irr_1.name_);
  micm::Index irr2_idx = state.variable_map_.at(irr_2.name_);

  std::vector<micm::Real> concentrations(5, 0);

  concentrations[a_idx] = 1.0;

  state.variables_[0] = concentrations;
  state.conditions_[0].temperature_ = 273;
  state.conditions_[0].pressure_ = 1000;

  for (micm::Index t{}; t < 100; ++t)
  {
    if (t > 50)
    {
      state.SetCustomRateParameter("r2", 0.0);
    }
    solver.UpdateStateParameters(state);
    std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
    micm::Real irr1 = state.variables_[0][irr1_idx];
    micm::Real irr2 = state.variables_[0][irr2_idx];
    auto result = solver.Solve(30.0, state);
    EXPECT_GE(state.variables_[0][irr1_idx], irr1);
    EXPECT_GE(state.variables_[0][irr2_idx], irr2);
  }
  std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
}
