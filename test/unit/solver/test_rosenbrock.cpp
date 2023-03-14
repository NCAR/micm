#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver.hpp>

#include <gtest/gtest.h>

TEST(ChapmanODESolver, DefaultConstructor){
  micm::RosenbrockSolver solver{};
}