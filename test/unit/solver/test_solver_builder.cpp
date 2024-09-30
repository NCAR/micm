#include <micm/solver/backward_euler.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#ifdef MICM_ENABLE_LLVM
  #include <micm/solver/jit/jit_solver_builder.hpp>
  #include <micm/solver/jit/jit_solver_parameters.hpp>
#endif
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

TEST(SolverBuilder, ThrowsMissingSystem)
{
  EXPECT_THROW(
      micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
          .SetNumberOfGridCells(1)
          .Build(),
      std::system_error);
}

TEST(SolverBuilder, ThrowsMissingReactions)
{
  EXPECT_THROW(
      micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
          .SetSystem(the_system)
          .SetNumberOfGridCells(1)
          .Build(),
      std::system_error);
  EXPECT_THROW(
      micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
          .SetSystem(the_system)
          .SetReactions({})
          .Build(),
      std::system_error);
}

TEST(SolverBuilder, CanBuildBackwardEuler)
{
  auto backward_euler = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
                            .SetSystem(the_system)
                            .SetReactions(reactions)
                            .SetNumberOfGridCells(1)
                            .Build();

  constexpr std::size_t L = 4;
  auto backward_euler_vector =
      micm::CpuSolverBuilder<
          micm::BackwardEulerSolverParameters,
          micm::VectorMatrix<double, L>,
          micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>(micm::BackwardEulerSolverParameters{})
          .SetSystem(the_system)
          .SetReactions(reactions)
          .SetNumberOfGridCells(1)
          .Build();
}

TEST(SolverBuilder, CanBuildRosenbrock)
{
  auto rosenbrock = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                        micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                        .SetSystem(the_system)
                        .SetReactions(reactions)
                        .SetNumberOfGridCells(1)
                        .Build();

  constexpr std::size_t L = 4;
  auto rosenbrock_vector = micm::CpuSolverBuilder<
                               micm::RosenbrockSolverParameters,
                               micm::VectorMatrix<double, L>,
                               micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>(
                               micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                               .SetSystem(the_system)
                               .SetReactions(reactions)
                               .SetNumberOfGridCells(1)
                               .Build();
}

/*
TEST(SolverBuilder, MismatchedToleranceSizeIsCaught)
{
  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  // too many
  params.absolute_tolerance_ = { 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
  EXPECT_ANY_THROW(micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(params)
                       .SetSystem(the_system)
                       .SetReactions(reactions)
                       .SetNumberOfGridCells(1)
                       .Build(););

  constexpr std::size_t L = 4;
  // too few
  params.absolute_tolerance_ = { 1e-6, 1e-6 };

  auto builder = micm::CpuSolverBuilder<
      micm::RosenbrockSolverParameters,
      micm::VectorMatrix<double, L>,
      micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>(params);

  EXPECT_ANY_THROW(builder.SetSystem(the_system).SetReactions(reactions).SetNumberOfGridCells(1).Build(););
}
*/