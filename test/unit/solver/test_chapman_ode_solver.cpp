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

TEST(ChapmanODESolver, simple_p_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 1);
  std::vector<double> number_densities(9, 1);
  double number_density_air{};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);

  for(auto& elem : forcing){
    std::cout << elem << std::endl;
  }

  // the truth values were calculated in fortran with old micm
  EXPECT_EQ(forcing[0], 0);
  EXPECT_EQ(forcing[1], 0);
  EXPECT_EQ(forcing[2], 0);
  EXPECT_EQ(forcing[3], 0);
  EXPECT_EQ(forcing[4], 0);
  EXPECT_EQ(forcing[5], -1);
  EXPECT_EQ(forcing[6], 3);
  EXPECT_EQ(forcing[7], 2);
  EXPECT_EQ(forcing[8], -2);
}

TEST(ChapmanODESolver, smaller_p_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 3e-8);
  std::vector<double> number_densities(9, 5e-6);
  double number_density_air{6e-14};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);

  // the truth values were calculated in fortran with old micm
  EXPECT_EQ(forcing[0], 0);
  EXPECT_EQ(forcing[1], 0);
  EXPECT_EQ(forcing[2], 0);
  EXPECT_EQ(forcing[3], 0);
  EXPECT_EQ(forcing[4], 0);
  EXPECT_NEAR(forcing[5], 1.49e-13, 0.01);
  EXPECT_NEAR(forcing[6], 4.55e-13, 0.01);
  EXPECT_NEAR(forcing[7], 1.5e-13, 0.01);
  EXPECT_NEAR(forcing[8], -3e-13, 0.01);
}