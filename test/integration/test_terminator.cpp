#include "terminator.hpp"

#include <micm/process/process.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

TEST(RosenbrockSolver, Terminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.relative_tolerance_ = 1.0e-8;
  parameters.max_number_of_steps_ = 100000;
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
}

template<std::size_t L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RosenbrockSolver, VectorTerminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.relative_tolerance_ = 1.0e-8;
  parameters.max_number_of_steps_ = 100000;
  {
    auto builder = VectorBuilder<1>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = VectorBuilder<2>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = VectorBuilder<3>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = VectorBuilder<4>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
}