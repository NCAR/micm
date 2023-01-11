#include <micm/solver/chapman_ode_solver_builder.hpp>
#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

TEST(SolverBuilder, DefaultConstructor){
  micm::ChapmanODESolverBuilder builder{};
}
