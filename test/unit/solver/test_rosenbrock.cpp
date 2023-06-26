#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver.hpp>
#include <micm/util/matrix.hpp>

TEST(ChapmanODESolver, DefaultConstructor)
{
  micm::RosenbrockSolver<micm::Matrix> solver{};
}