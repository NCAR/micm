#include <micm/solver/chapman_ode_solver_builder.hpp>
#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

TEST(SolverBuilder, DefaultConstructor){
  micm::ChapmanODESolverBuilder builder{};
}

TEST(SolverBuilder, ForProcess){
  micm::ChapmanODESolverBuilder builder{};
  builder.For(micm::Process());
}

TEST(SolverBuilder, VectorProcess){
  micm::ChapmanODESolverBuilder builder{};
  builder.For(
    std::vector<micm::Process>{micm::Process(), micm::Process()}
  );
}

TEST(SolverBuilder, Build){
  micm::ChapmanODESolverBuilder builder{};
  builder.Build();
}
