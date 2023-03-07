#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/solver.hpp>
#include <algorithm>
#include <random>

#include <gtest/gtest.h>

static const double absolute_tolerance = 1e-4;

void TestDefaultConstructor(micm::RosenbrockSolver& solver){
  EXPECT_EQ(solver.parameters_.stages_, 3);
}
TEST(ChapmanMechanismHardCodedAndGeneral, DefaultConstructor){
  micm::ChapmanODESolver hard_coded{};
  // micm::RosenbrockSolver general{};
  TestDefaultConstructor(hard_coded);
  // TestDefaultConstructor(general);
}

void TestLinSolve(micm::RosenbrockSolver& solver){
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
TEST(ChapmanMechanismHardCodedAndGeneral, lin_solve){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestLinSolve(solver);
  // TestLinSolve(general);
}

void TestReactionNames(micm::RosenbrockSolver& solver){
  auto names = solver.reaction_names();
  ASSERT_EQ(names.size(), 7);
}
TEST(ChapmanMechanismHardCodedAndGeneral, ReactionNames){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestReactionNames(solver);
}

void TestPhotolysisNames(micm::RosenbrockSolver& solver){
  auto names = solver.photolysis_names();
  ASSERT_EQ(names.size(), 3);
}
TEST(ChapmanMechanismHardCodedAndGeneral, PhotolysisNames){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestPhotolysisNames(solver);
}

void TestSpeciesNames(micm::RosenbrockSolver& solver) {
  auto names = solver.species_names();
  ASSERT_EQ(names.size(), 9);
}
TEST(ChapmanMechanismHardCodedAndGeneral, SpeciesNames){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSpeciesNames(solver);
}

void TestSimpleForce(micm::RosenbrockSolver& solver){
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
TEST(ChapmanMechanismHardCodedAndGeneral, simple_force){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSimpleForce(solver);
}

void TestSmallerForce(micm::RosenbrockSolver& solver){
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
TEST(ChapmanMechanismHardCodedAndGeneral, smaller_force){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSmallerForce(solver);
}

void TestFactoredAlphaMinusJac(micm::RosenbrockSolver& solver){
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
TEST(ChapmanMechanismHardCodedAndGeneral, factored_alpha_minus_jac){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestFactoredAlphaMinusJac(solver);
}

void TestDforceDyTimesVector(micm::RosenbrockSolver& solver){
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
TEST(ChapmanMechanismHardCodedAndGeneral, dforce_dy_times_vector){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestDforceDyTimesVector(solver);
}

void TestSolve(micm::RosenbrockSolver& solver){
  std::vector<double> number_densities = { 1,    3.92e-1, 1.69e-2, 0,     3.29e1, 0,     0,   8.84, 0};
                                         //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  solver.calculate_rate_constants(temperature, pressure);

  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);
  EXPECT_EQ(results.state_, micm::Solver::SolverState::Converged);
  EXPECT_NEAR(results.result_[0], 1, absolute_tolerance);
  EXPECT_NEAR(results.result_[1], 0.392, absolute_tolerance);
  EXPECT_NEAR(results.result_[2], 0.0169, absolute_tolerance);
  EXPECT_NEAR(results.result_[3], 0, absolute_tolerance);
  EXPECT_NEAR(results.result_[4], 32.9, absolute_tolerance);
  EXPECT_NEAR(results.result_[5], 1.8039e-41, absolute_tolerance);
  EXPECT_NEAR(results.result_[6], 0.00176789, absolute_tolerance);
  EXPECT_NEAR(results.result_[7], 8.83912, absolute_tolerance);
  EXPECT_NEAR(results.result_[8], 4.5031e-36, absolute_tolerance);
}
TEST(ChapmanMechanismHardCodedAndGeneral, Solve){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSolve(solver);
}

void TestSolve10TimesLarger(micm::RosenbrockSolver& solver){
  std::vector<double> number_densities = { 1,    3.92e-1, 1.69e-2, 0,     3.29e1, 0,     0,   8.84, 0};
                                         //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  for(auto& elem: number_densities){
    elem *= 10;
  }

  solver.calculate_rate_constants(temperature, pressure);
  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);
  EXPECT_NEAR(results.result_[0], 10, absolute_tolerance);
  EXPECT_NEAR(results.result_[1], 3.92, absolute_tolerance);
  EXPECT_NEAR(results.result_[2], 0.169, absolute_tolerance);
  EXPECT_NEAR(results.result_[3], 0, absolute_tolerance);
  EXPECT_NEAR(results.result_[4], 329, absolute_tolerance);
  EXPECT_NEAR(results.result_[5], 1.8039e-38, absolute_tolerance);
  EXPECT_NEAR(results.result_[6], 0.0176789, absolute_tolerance);
  EXPECT_NEAR(results.result_[7], 88.3912, absolute_tolerance);
  EXPECT_NEAR(results.result_[8], 4.5031e-33, absolute_tolerance);
}
TEST(ChapmanMechanismHardCodedAndGeneral, solve_10_times_larger){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSolve10TimesLarger(solver);
}

void TestSolve10TimesSmaller(micm::RosenbrockSolver& solver){
  std::vector<double> number_densities = { 1,    3.92e-1, 1.69e-2, 0,     3.29e1, 0,     0,   8.84, 0};
                                         //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  for(auto& elem: number_densities){
    elem /= 10;
  }

  solver.calculate_rate_constants(temperature, pressure);
  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);
  EXPECT_NEAR(results.result_[0], 0.1, absolute_tolerance);
  EXPECT_NEAR(results.result_[1], 0.0392, absolute_tolerance);
  EXPECT_NEAR(results.result_[2], 0.00169, absolute_tolerance);
  EXPECT_NEAR(results.result_[3], 0, absolute_tolerance);
  EXPECT_NEAR(results.result_[4], 3.29, absolute_tolerance);
  EXPECT_NEAR(results.result_[5], 1.8039e-44, absolute_tolerance);
  EXPECT_NEAR(results.result_[6], 0.000176789, absolute_tolerance);
  EXPECT_NEAR(results.result_[7], 0.883912, absolute_tolerance);
  EXPECT_NEAR(results.result_[8], 4.5031e-39, absolute_tolerance);
}
TEST(RegressionChapmanODESolver, solve_10_times_smaller){
  micm::ChapmanODESolver solver{};
  // micm::RosenbrockSolver general{};
  TestSolve10TimesSmaller(solver);
}

void TestSolveWithRandomNumberDensities(micm::RosenbrockSolver& solver){
  std::vector<double> number_densities(9);
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  std::generate(number_densities.begin(), number_densities.end(),
                [&] { return distribution(generator); });

  solver.calculate_rate_constants(temperature, pressure);

  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);
  EXPECT_NEAR(results.result_[0], 7.8259e-06, absolute_tolerance);
  EXPECT_NEAR(results.result_[1], 0.131538, absolute_tolerance);
  EXPECT_NEAR(results.result_[2], 0.755605, absolute_tolerance);
  EXPECT_NEAR(results.result_[3], 0.45865, absolute_tolerance);
  EXPECT_NEAR(results.result_[4], 0.532767, absolute_tolerance);
  EXPECT_NEAR(results.result_[5], 0.218966, absolute_tolerance);
  EXPECT_NEAR(results.result_[6], 0.0471811, absolute_tolerance);
  EXPECT_NEAR(results.result_[7], 0.678804, absolute_tolerance);
  EXPECT_NEAR(results.result_[8], 0.679289, absolute_tolerance);
}
TEST(ChapmanMechanismHardCodedAndGeneral, solve_with_random_number_densities){
  micm::ChapmanODESolver solver{};
  TestSolveWithRandomNumberDensities(solver);
}