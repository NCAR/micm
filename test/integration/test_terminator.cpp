#include "terminator.hpp"

#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(RosenbrockSolver, Terminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
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

template<micm::Index L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<micm::Real, L>,
    micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RosenbrockSolver, VectorTerminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
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