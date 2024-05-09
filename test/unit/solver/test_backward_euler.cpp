#include <micm/solver/backward_euler.hpp>
#include <micm/solver/rosenbrock.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

TEST(BackwardEuler, DefaultConstructor)
{
  EXPECT_NO_THROW(micm::BackwardEuler be);
}

TEST(BackwardEuler, CanCallSolve)
{
  micm::BackwardEuler be;
  double time_step = 1.0;

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yields(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yields(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  options.ignore_unused_species_ = true;

  micm::RosenbrockSolver<> rosenbrock{ micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
                                       std::vector<micm::Process>{ r1, r2 },
                                       options };

  auto state = rosenbrock.GetState();

  state.variables_[0] = { 1.0, 0.0, 0.0 };
  state.conditions_[0].temperature_ = 272.5;
  state.conditions_[0].pressure_ = 101253.3;
  state.conditions_[0].air_density_ = 1e6;

  auto linear_solver = rosenbrock.linear_solver_;
  auto process_set = rosenbrock.process_set_;
  auto processes = std::vector<micm::Process>{ r1, r2 };

  EXPECT_NO_THROW(be.Solve(
      time_step, state, linear_solver, process_set, processes, rosenbrock.state_parameters_.jacobian_diagonal_elements_));
}
