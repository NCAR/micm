#include <micm/GPU.hpp>
#include <micm/Process.hpp>
#include <micm/Solver.hpp>
#include <micm/Util.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

template<std::size_t L>
using GpuBuilder = micm::CudaSolverBuilderInPlace<micm::CudaRosenbrockSolverParameters, L>;

namespace
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase)
                         .Build();
  micm::System the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(SolverBuilder, CanBuildCudaSolvers)
{
  constexpr std::size_t L = 4;
  auto cuda_rosenbrock = GpuBuilder<L>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                             .SetSystem(the_system)
                             .SetReactions(reactions)
                             .Build();
  auto state = cuda_rosenbrock.GetState(L);
}