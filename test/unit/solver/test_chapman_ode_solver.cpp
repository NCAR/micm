#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

TEST(SolverBuilder, DefaultConstructor){
  micm::ChapmanODESolver solver{};
}
