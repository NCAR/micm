#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

using SparseMatrixTest = micm::SparseMatrix<double>;

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto a = micm::Species("a");
  auto b = micm::Species("b");
  auto c = micm::Species("c");
  auto irr_1 = micm::Species("irr_1");
  auto irr_2 = micm::Species("irr_2");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c, irr_1, irr_2 } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ Yields(b, 1), Yields(irr_1, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 =
      micm::Process::Create()
          .SetReactants({ b })
          .SetProducts({ Yields(c, 1), Yields(irr_2, 1) })
          .SetRateConstant(micm::UserDefinedRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" }))
          .SetPhase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::vector<micm::Process>{ r1, r2 }, options
  };

  auto state = solver.GetState();

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
    std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
    auto result = solver.Solve(30.0, state);
    EXPECT_GE(result.result_[0][irr1_idx], state.variables_[0][irr1_idx]);
    EXPECT_GE(result.result_[0][irr2_idx], state.variables_[0][irr2_idx]);
    state.variables_ = result.result_;
  }
  std::cout << state.variables_[0][irr1_idx] << " " << state.variables_[0][irr2_idx] << std::endl;
}
