#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/solver.hpp>

#include <gtest/gtest.h>

TEST(ChapmanODESolver, DefaultConstructor){
  micm::ChapmanODESolver solver{};

  EXPECT_EQ(solver.parameters_.stages_, 3);
}

TEST(ChapmanODESolver, lin_solve){
  micm::ChapmanODESolver solver{};
  std::vector<double> jacobian(23, 1), b(9, 0.5);
  auto solved = solver.lin_solve(b, jacobian);

  EXPECT_EQ(solved[0], 0.5);
  EXPECT_EQ(solved[1], 0.5);
  EXPECT_EQ(solved[2], 0.5);
  EXPECT_EQ(solved[3], 0.5);
  EXPECT_EQ(solved[4], 0.5);
  EXPECT_EQ(solved[5], -.5);
  EXPECT_EQ(solved[6], -1);
  EXPECT_EQ(solved[7], 0.5);
  EXPECT_EQ(solved[8], 0);
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

TEST(ChapmanODESolver, simple_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 1);
  std::vector<double> number_densities(9, 1);
  double number_density_air{};

  auto forcing = solver.force(rate_constants, number_densities, number_density_air);

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

TEST(ChapmanODESolver, smaller_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 3e-8);
  std::vector<double> number_densities(9, 5e-6);
  double number_density_air{6e-14};

  auto forcing = solver.force(rate_constants, number_densities, number_density_air);

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

  auto jacobian = solver.factored_alpha_minus_jac(dforce_dy, alpha);

  // the truth values were calculated in fortran with old micm
  EXPECT_NEAR(jacobian[0], 1.000, 0.01);
  EXPECT_NEAR(jacobian[1], -1.000, 0.01);
  EXPECT_NEAR(jacobian[2], -1.000, 0.01);
  EXPECT_NEAR(jacobian[3], -1.000, 0.01);
  EXPECT_NEAR(jacobian[4], 1.000, 0.01);
  EXPECT_NEAR(jacobian[5], 1.000, 0.01);
  EXPECT_NEAR(jacobian[6], 1.000, 0.01);
  EXPECT_NEAR(jacobian[7], 1.000, 0.01);
  EXPECT_NEAR(jacobian[8], -1.000, 0.01);
  EXPECT_NEAR(jacobian[9], -1.000, 0.01);
  EXPECT_NEAR(jacobian[10], 1.000, 0.01);
  EXPECT_NEAR(jacobian[11], -1.000, 0.01);
  EXPECT_NEAR(jacobian[12], 1.000, 0.01);
  EXPECT_NEAR(jacobian[13], -1.000, 0.01);
  EXPECT_NEAR(jacobian[14], -1.000, 0.01);
  EXPECT_NEAR(jacobian[15], -1.000, 0.01);
  EXPECT_NEAR(jacobian[16], -2.000, 0.01);
  EXPECT_NEAR(jacobian[17], -1.000, 0.01);
  EXPECT_NEAR(jacobian[18], 3.000, 0.01);
  EXPECT_NEAR(jacobian[19], -1.000, 0.01);
  EXPECT_NEAR(jacobian[20], -2.000, 0.01);
  EXPECT_NEAR(jacobian[21], -3.000, 0.01);
  EXPECT_NEAR(jacobian[22], 0.125, 0.01);
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

TEST(ChapmanODESolver, Solve){
  micm::ChapmanODESolver solver{};
  std::vector<double> number_densities = { 1,    3.92e-1, 1.69e-2, 0,     3.29e1, 0,     0,   8.84, 0};
                                         //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  solver.calculate_rate_constants(temperature, pressure);

  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);

  std::cout << "solver state: " << micm::state_to_string(results.state_) << "\n";
}