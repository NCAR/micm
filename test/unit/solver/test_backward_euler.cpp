#include <micm/solver/solver_builder.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

namespace
{
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

  auto the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(BackwardEuler, CanCallSolve)
{
  auto params = micm::BackwardEulerSolverParameters();
  params.absolute_tolerance_ = { 1e-6, 1e-6, 1e-6 };

  auto be = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(params)
                .SetSystem(the_system)
                .SetReactions(reactions)
                .SetNumberOfGridCells(1)
                .Build();
  double time_step = 1.0;

  auto state = be.GetState();

  state.variables_[0] = { 1.0, 0.0, 0.0 };
  state.conditions_[0].temperature_ = 272.5;
  state.conditions_[0].pressure_ = 101253.3;
  state.conditions_[0].air_density_ = 1e6;

  EXPECT_NO_THROW(auto result = be.Solve(time_step, state));
}
