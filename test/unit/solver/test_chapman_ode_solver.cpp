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

TEST(ChapmanODESolver, factored_alpha_minus_jac){
  micm::ChapmanODESolver solver{};
  std::vector<double> dforce_dy(23, 1);
  double alpha{2};

  auto LU = solver.factored_alpha_minus_jac(dforce_dy, alpha);

  // the truth values were calculated in fortran with old micm
  EXPECT_NEAR(LU[0], 1.000, 0.01);
  EXPECT_NEAR(LU[1], -1.000, 0.01);
  EXPECT_NEAR(LU[2], -1.000, 0.01);
  EXPECT_NEAR(LU[3], -1.000, 0.01);
  EXPECT_NEAR(LU[4], 1.000, 0.01);
  EXPECT_NEAR(LU[5], 1.000, 0.01);
  EXPECT_NEAR(LU[6], 1.000, 0.01);
  EXPECT_NEAR(LU[7], 1.000, 0.01);
  EXPECT_NEAR(LU[8], -1.000, 0.01);
  EXPECT_NEAR(LU[9], -1.000, 0.01);
  EXPECT_NEAR(LU[10], 1.000, 0.01);
  EXPECT_NEAR(LU[11], -1.000, 0.01);
  EXPECT_NEAR(LU[12], 1.000, 0.01);
  EXPECT_NEAR(LU[13], -1.000, 0.01);
  EXPECT_NEAR(LU[14], -1.000, 0.01);
  EXPECT_NEAR(LU[15], -1.000, 0.01);
  EXPECT_NEAR(LU[16], -2.000, 0.01);
  EXPECT_NEAR(LU[17], -1.000, 0.01);
  EXPECT_NEAR(LU[18], 3.000, 0.01);
  EXPECT_NEAR(LU[19], -1.000, 0.01);
  EXPECT_NEAR(LU[20], -2.000, 0.01);
  EXPECT_NEAR(LU[21], -3.000, 0.01);
  EXPECT_NEAR(LU[22], 0.125, 0.01);
}

TEST(ChapmanODESolver, dforce_dy_time_vector){
  micm::ChapmanODESolver solver{};
  std::vector<double> dforce_dy(23, 1);
  std::vector<double> vector(23, 0.5);

  auto product = solver.dforce_dy_times_vector(dforce_dy, vector);

  // the truth values were calculated in fortran with old micm
  EXPECT_NEAR(product[0], 0, 0.01);
  EXPECT_NEAR(product[1], 0, 0.01);
  EXPECT_NEAR(product[2], 0, 0.01);
  EXPECT_NEAR(product[3], 0, 0.01);
  EXPECT_NEAR(product[4], 0, 0.01);
  EXPECT_NEAR(product[5], 2, 0.01);
  EXPECT_NEAR(product[6], 3, 0.01);
  EXPECT_NEAR(product[7], 2, 0.01);
  EXPECT_NEAR(product[8], 2, 0.01);
  EXPECT_NEAR(product[9], 0, 0.01);
  EXPECT_NEAR(product[10], 0, 0.01);
  EXPECT_NEAR(product[11], 0, 0.01);
  EXPECT_NEAR(product[12], 0, 0.01);
  EXPECT_NEAR(product[13], 0, 0.01);
  EXPECT_NEAR(product[14], 0, 0.01);
  EXPECT_NEAR(product[15], 0, 0.01);
  EXPECT_NEAR(product[16], 0, 0.01);
  EXPECT_NEAR(product[17], 0, 0.01);
  EXPECT_NEAR(product[18], 0, 0.01);
  EXPECT_NEAR(product[19], 0, 0.01);
  EXPECT_NEAR(product[20], 0, 0.01);
  EXPECT_NEAR(product[21], 0, 0.01);
  EXPECT_NEAR(product[22], 0, 0.01);
}