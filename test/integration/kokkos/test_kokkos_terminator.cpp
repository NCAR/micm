#include "../terminator.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(RosenbrockSolver, Terminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.max_number_of_steps_ = 100000;
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
}

template<micm::Index L>
using VectorBuilder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters, L>;

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

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
