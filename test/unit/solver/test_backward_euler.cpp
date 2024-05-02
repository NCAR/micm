#include <micm/solver/backward_euler.hpp>
#include <micm/solver/rosenbrock.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

TEST(BackwardEuler, DefaultConstructor)
{
  EXPECT_NO_THROW(micm::BackwardEuler be);
}

TEST(BackwardEuler, CanCallSolve) {
  micm::BackwardEuler be;
  double time_step = 1.0;

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(b, 1) })
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(c, 1) })
                         .rate_constant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .phase(gas_phase);

  auto options = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters();
  options.ignore_unused_species_ = true;

  micm::RosenbrockSolver<> rosenbrock{
    micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<micm::Process>{ r1, r2 },
    options
  };

  auto state = rosenbrock.GetState();
  auto linear_solver = rosenbrock.linear_solver_;
  auto process_set = rosenbrock.process_set_;
  auto processes = std::vector<micm::Process>{ r1, r2 };

  EXPECT_NO_THROW(be.Solve(time_step, state, linear_solver, process_set, processes, rosenbrock.state_parameters_.jacobian_diagonal_elements_));
}
