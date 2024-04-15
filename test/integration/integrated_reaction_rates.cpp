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

TEST(ChapmanIntegration, CanBuildChapmanSystem)
{
  auto a = micm::Species("a");
  auto b = micm::Species("b");
  auto c = micm::Species("c");
  auto irr_1  = micm::Species("irr_1");
  auto irr_2 = micm::Species("irr_2");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c, irr_1, irr_2 } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1), yields(irr_1, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1), yields(irr_2, 1) })
                         .rate_constant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .phase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters();

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest> solver{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2  },
    options
  };

  auto state = solver.GetState();

  std::vector<double> concentrations{ 1, 0, 0, 0, 0 };
  state.variables_[0] = concentrations;
  state.conditions_[0].temperature_ = 273;
  state.conditions_[0].pressure_ = 1000;

  for (double t{}; t < 100; ++t)
  {
    std::cout << state.variables_[0][state.variable_map_.at(irr_1.name_)] << " " << state.variables_[0][state.variable_map_.at(irr_2.name_)] << std::endl;
    auto result = solver.Solve(30.0, state);
    EXPECT_GE(result.result_[0][state.variable_map_.at(irr_1.name_)], state.variables_[0][state.variable_map_.at(irr_1.name_)]);
    EXPECT_GE(result.result_[0][state.variable_map_.at(irr_1.name_)], state.variables_[0][state.variable_map_.at(irr_1.name_)]);
    state.variables_ = result.result_;
  }
  std::cout << state.variables_[0][state.variable_map_.at(irr_1.name_)] << " " << state.variables_[0][state.variable_map_.at(irr_2.name_)] << std::endl;
}
