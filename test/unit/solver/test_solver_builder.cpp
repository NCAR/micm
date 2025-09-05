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
                         .SetProducts({ micm::Yield(b, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .SetPhase(gas_phase);

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ micm::Yield(c, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstant(
                             micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .SetPhase(gas_phase);
  micm::System the_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  std::vector<micm::Process> reactions = { r1, r2 };
}  // namespace

TEST(SolverBuilder, ThrowsMissingSystem)
{
  EXPECT_THROW(
      micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{}).Build(),
      std::system_error);
}

TEST(SolverBuilder, ThrowsMissingReactions)
{
  EXPECT_THROW(
      micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
          .SetSystem(the_system)
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
                            .Build();

  constexpr std::size_t L = 4;
  auto backward_euler_vector =
      micm::CpuSolverBuilder<
          micm::BackwardEulerSolverParameters,
          micm::VectorMatrix<double, L>,
          micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>(micm::BackwardEulerSolverParameters{})
          .SetSystem(the_system)
          .SetReactions(reactions)
          .Build();
  EXPECT_EQ(backward_euler_vector.GetSystem().gas_phase_.name_, the_system.gas_phase_.name_);
  EXPECT_EQ(backward_euler_vector.GetSystem().gas_phase_.species_.size(), the_system.gas_phase_.species_.size());
  EXPECT_EQ(backward_euler_vector.GetSystem().phases_.size(), the_system.phases_.size());
  EXPECT_GT(backward_euler.MaximumNumberOfGridCells(), 1e8);
  EXPECT_EQ(backward_euler_vector.MaximumNumberOfGridCells(), 4);
}

TEST(SolverBuilder, CanBuildRosenbrock)
{
  auto rosenbrock = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                        micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                        .SetSystem(the_system)
                        .SetReactions(reactions)
                        .Build();

  constexpr std::size_t L = 4;
  auto rosenbrock_vector = micm::CpuSolverBuilder<
                               micm::RosenbrockSolverParameters,
                               micm::VectorMatrix<double, L>,
                               micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>(
                               micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                               .SetSystem(the_system)
                               .SetReactions(reactions)
                               .Build();

  EXPECT_EQ(rosenbrock_vector.GetSystem().gas_phase_.name_, the_system.gas_phase_.name_);
  EXPECT_EQ(rosenbrock_vector.GetSystem().gas_phase_.species_.size(), the_system.gas_phase_.species_.size());
  EXPECT_EQ(rosenbrock_vector.GetSystem().phases_.size(), the_system.phases_.size());
  EXPECT_GT(rosenbrock.MaximumNumberOfGridCells(), 1e8);
  EXPECT_EQ(rosenbrock_vector.MaximumNumberOfGridCells(), 4);
}

TEST(SolverBuilder, CanBuildBackwardEulerOverloadedSolverMethod)
{
  auto solver = micm::CpuSolverBuilder<micm::BackwardEulerSolverParameters>(micm::BackwardEulerSolverParameters{})
                    .SetSystem(the_system)
                    .SetReactions(reactions)
                    .Build();
  auto state = solver.GetState(1);
  auto options = micm::BackwardEulerSolverParameters();
  auto solve = solver.Solve(5, state, options);

  EXPECT_EQ(solve.final_time_, 5);
  EXPECT_EQ(solve.stats_.function_calls_, 2);
  EXPECT_EQ(solve.stats_.jacobian_updates_, 2);
  EXPECT_EQ(solve.stats_.number_of_steps_, 2);
  EXPECT_EQ(solve.stats_.solves_, 2);

  options.small_ = 1.0;
  options.max_number_of_steps_ = 1.0;

  solve = solver.Solve(5, state, options);

  EXPECT_EQ(solve.final_time_, 0.03125);
  EXPECT_EQ(solve.stats_.function_calls_, 6);
  EXPECT_EQ(solve.stats_.jacobian_updates_, 6);
  EXPECT_EQ(solve.stats_.number_of_steps_, 6);
  EXPECT_EQ(solve.stats_.solves_, 6);
}

TEST(SolverBuilder, CanBuildRosenbrockOverloadedSolveMethod)
{
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                    .SetSystem(the_system)
                    .SetReactions(reactions)
                    .Build();
  auto state = solver.GetState(1);
  state.variables_[0] = { 1.0, 0.0, 0.0 };

  auto solve = solver.Solve(5, state);

  EXPECT_EQ(solve.final_time_, 5);
  EXPECT_EQ(solve.stats_.function_calls_, 18);
  EXPECT_EQ(solve.stats_.jacobian_updates_, 9);
  EXPECT_EQ(solve.stats_.number_of_steps_, 9);
  EXPECT_EQ(solve.stats_.solves_, 27);

  options.h_min_ = 15.0;
  options.max_number_of_steps_ = 6.0;

  state.variables_[0] = { 1.0, 0.0, 0.0 };
  solve = solver.Solve(5, state, options);

  EXPECT_EQ(solve.final_time_, 5);
  EXPECT_EQ(solve.stats_.function_calls_, 2);
  EXPECT_EQ(solve.stats_.jacobian_updates_, 1);
  EXPECT_EQ(solve.stats_.number_of_steps_, 1);
  EXPECT_EQ(solve.stats_.solves_, 3);

  EXPECT_EQ(solver.solver_parameters_.h_min_, 15.0);
  EXPECT_EQ(solver.solver_parameters_.max_number_of_steps_, 6.0);
}