#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver.hpp>

TEST(ChapmanODESolver, DefaultConstructor)
{
  micm::RosenbrockSolver solver{};
}