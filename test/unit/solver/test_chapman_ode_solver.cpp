#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

TEST(SolverBuilder, DefaultConstructor){
  micm::ChapmanODESolver solver{};
}

TEST(SolverBuilder, Solve){
  micm::ChapmanODESolver solver{};
  double state[] = {1, 2, 3, 4};
  solver.Solve(state);
}

