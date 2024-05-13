#include <micm/solver/backward_euler.hpp>
#include <micm/solver/solver_builder.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

namespace {
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
  micm::System the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}

TEST(SolverBuilder, CanBuildBackwardEuler)
{
  auto backward_euler = micm::SolverBuilder()
                            .SetSystem(the_system)
                            .SetReactions(reactions)
                            .SetNumberOfGridCells(1)
                            .BackwardEuler(micm::BackwardEulerSolverParameters{})
                            .Build();
}