#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

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
  micm::System the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(SolverBuilder, CanBuildJitRosenbrock)
{
  constexpr std::size_t L = 4;
  auto jit_rosenbrock = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>(
                            micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                            .SetSystem(the_system)
                            .SetReactions(reactions)
                            .SetNumberOfGridCells(L)
                            .Build();
}

TEST(SolverBuilder, MismatchedToleranceSizeIsCaught)
{
  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // too many
  params.absolute_tolerance_ = { 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };

  constexpr std::size_t L = 4;
  using jit_builder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>;

  EXPECT_ANY_THROW(jit_builder(params).SetSystem(the_system).SetReactions(reactions).SetNumberOfGridCells(L).Build(););

  // too few
  params.absolute_tolerance_ = { 1e-6, 1e-6 };
  EXPECT_ANY_THROW(jit_builder(params).SetSystem(the_system).SetReactions(reactions).SetNumberOfGridCells(L).Build(););
}