#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

TEST(ChapmanODESolver, DefaultConstructor){
  micm::ChapmanODESolver solver{};
}

TEST(ChapmanODESolver, Solve){
  micm::ChapmanODESolver solver{};
  double state[] = {1, 2, 3, 4};
  solver.Solve(state);
}

TEST(ChapmanODESolver, ReactionNames){
  micm::ChapmanODESolver solver{};
  auto names = solver.reaction_names();
  ASSERT_EQ(names.size(), 7);
}

TEST(ChapmanODESolver, PhotolysisNames){
  micm::ChapmanODESolver solver{};
  auto names = solver.photolysis_names();
  ASSERT_EQ(names.size(), 3);
}

TEST(ChapmanODESolver, SpeciesNames){
  micm::ChapmanODESolver solver{};
  auto names = solver.species_names();
  ASSERT_EQ(names.size(), 9);
}

TEST(ChapmanODESolver, p_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants{};
  std::vector<double> number_densities{};
  double number_density_air{};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);
}