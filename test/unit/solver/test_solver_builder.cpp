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

TEST(SolverBuilder, CanBuildBackwardEulerOverloadedSolverMethod)
{
  auto solver = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
                            .SetSystem(the_system)
                            .SetReactions(reactions)
                            .SetNumberOfGridCells(1)
                            .Build();
  auto state = solver.GetState();
  auto options = micm::BackwardEulerSolverParameters();
  auto solve = solver.Solve(5, state, options);

  ASSERT_EQ(solve.final_time_, 5);
  ASSERT_EQ(solve.stats_.function_calls_, 2);
  ASSERT_EQ(solve.stats_.jacobian_updates_, 2);
  ASSERT_EQ(solve.stats_.number_of_steps_, 2);
  ASSERT_EQ(solve.stats_.solves_, 2);
  
  options.small_ = 1.0;
  options.max_number_of_steps_ = 1.0;

  solve = solver.Solve(5, state, options);

  ASSERT_EQ(solve.final_time_, 0.03125);
  ASSERT_EQ(solve.stats_.function_calls_, 6);
  ASSERT_EQ(solve.stats_.jacobian_updates_, 6);
  ASSERT_EQ(solve.stats_.number_of_steps_, 6);
  ASSERT_EQ(solve.stats_.solves_, 6);
}

TEST(SolverBuilder, CanBuildRosenbrockOverloadedSolveMethod)
{
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                        .SetSystem(the_system)
                        .SetReactions(reactions)
                        .SetNumberOfGridCells(1)
                        .Build();
  auto state = solver.GetState();
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solve = solver.Solve(5, state, options);

  ASSERT_EQ(solve.final_time_, 5);
  ASSERT_EQ(solve.stats_.function_calls_, 20);
  ASSERT_EQ(solve.stats_.jacobian_updates_, 10);
  ASSERT_EQ(solve.stats_.number_of_steps_, 10);
  ASSERT_EQ(solve.stats_.solves_, 30);
  
  options.h_min_ = 15.0;
  options.max_number_of_steps_ = 6.0;

  solve = solver.Solve(5, state, options);

  ASSERT_EQ(solve.final_time_, 5);
  ASSERT_EQ(solve.stats_.function_calls_, 2);
  ASSERT_EQ(solve.stats_.jacobian_updates_, 1);
  ASSERT_EQ(solve.stats_.number_of_steps_, 1);
  ASSERT_EQ(solve.stats_.solves_, 3);
}