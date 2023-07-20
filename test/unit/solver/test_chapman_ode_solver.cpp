#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "../../regression/RosenbrockChapman/chapman_ode_solver.hpp"

static const double absolute_tolerance = 1e-4;

void TestDefaultConstructor(micm::ChapmanODESolver& solver)
{
  EXPECT_EQ(solver.parameters_.stages_, 3);
}
TEST(ChapmanMechanismHardCodedAndGeneral, DefaultConstructor)
{
  micm::ChapmanODESolver hard_coded{};
  // micm::ChapmanODESolver general{};
  TestDefaultConstructor(hard_coded);
  // TestDefaultConstructor(general);
}

void TestLinSolve(micm::ChapmanODESolver& solver)
{
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
TEST(ChapmanMechanismHardCodedAndGeneral, lin_solve)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestLinSolve(solver);
  // TestLinSolve(general);
}

void TestReactionNames(micm::ChapmanODESolver& solver)
{
  auto names = solver.reaction_names();
  ASSERT_EQ(names.size(), 7);
}
TEST(ChapmanMechanismHardCodedAndGeneral, ReactionNames)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestReactionNames(solver);
}

void TestPhotolysisNames(micm::ChapmanODESolver& solver)
{
  auto names = solver.photolysis_names();
  ASSERT_EQ(names.size(), 3);
}
TEST(ChapmanMechanismHardCodedAndGeneral, PhotolysisNames)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestPhotolysisNames(solver);
}

void TestSpeciesNames(micm::ChapmanODESolver& solver)
{
  auto names = solver.species_names();
  ASSERT_EQ(names.size(), 9);
}
TEST(ChapmanMechanismHardCodedAndGeneral, SpeciesNames)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSpeciesNames(solver);
}

void TestSimpleForce(micm::ChapmanODESolver& solver)
{
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
TEST(ChapmanMechanismHardCodedAndGeneral, simple_force)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSimpleForce(solver);
}

void TestSmallerForce(micm::ChapmanODESolver& solver)
{
  std::vector<double> rate_constants(9, 3e-8);
  std::vector<double> number_densities(9, 5e-6);
  double number_density_air{ 6e-14 };

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
TEST(ChapmanMechanismHardCodedAndGeneral, smaller_force)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSmallerForce(solver);
}

void TestFactoredAlphaMinusJac(micm::ChapmanODESolver& solver)
{
  std::vector<double> dforce_dy(23, 1);
  double alpha{ 2 };

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
TEST(ChapmanMechanismHardCodedAndGeneral, factored_alpha_minus_jac)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestFactoredAlphaMinusJac(solver);
}

void TestDforceDyTimesVector(micm::ChapmanODESolver& solver)
{
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
TEST(ChapmanMechanismHardCodedAndGeneral, dforce_dy_times_vector)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestDforceDyTimesVector(solver);
}

void TestSolve(micm::ChapmanODESolver& solver)
{
  micm::State<micm::Matrix> state = solver.GetState();
  state.variables_[0] = { 1, 3.92e-1, 1.69e-2, 0, 3.29e1, 0, 0, 8.84, 0 };
  //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  state.conditions_[0].temperature_ = 273.15;
  state.conditions_[0].pressure_ = 1000 * 100;  // 1000 hPa
  state.conditions_[0].air_density_ = 2.7e19;
  state.custom_rate_parameters_[0][0] = 1.0e-4;
  state.custom_rate_parameters_[0][1] = 1.0e-5;
  state.custom_rate_parameters_[0][2] = 1.0e-6;
  double time_start = 0;
  double time_end = 1;

  solver.UpdateState(state);

  auto results = solver.Solve(time_start, time_end, state);
  EXPECT_EQ(results.state_, micm::ChapmanODESolver::SolverState::Converged);
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
TEST(ChapmanMechanismHardCodedAndGeneral, Solve)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSolve(solver);
}

void TestSolve10TimesLarger(micm::ChapmanODESolver& solver)
{
  micm::State<micm::Matrix> state = solver.GetState();
  state.variables_[0] = { 1, 3.92e-1, 1.69e-2, 0, 3.29e1, 0, 0, 8.84, 0 };
  //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  state.conditions_[0].temperature_ = 273.15;
  state.conditions_[0].pressure_ = 1000 * 100;  // 1000 hPa
  state.conditions_[0].air_density_ = 2.7e19;
  state.custom_rate_parameters_[0][0] = 1.0e-4;
  state.custom_rate_parameters_[0][1] = 1.0e-5;
  state.custom_rate_parameters_[0][2] = 1.0e-6;
  double time_start = 0;
  double time_end = 1;

  for (auto& elem : state.variables_[0])
  {
    elem *= 10;
  }

  solver.UpdateState(state);
  auto results = solver.Solve(time_start, time_end, state);
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
TEST(ChapmanMechanismHardCodedAndGeneral, solve_10_times_larger)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSolve10TimesLarger(solver);
}

void TestSolve10TimesSmaller(micm::ChapmanODESolver& solver)
{
  micm::State<micm::Matrix> state = solver.GetState();
  state.variables_[0] = { 1, 3.92e-1, 1.69e-2, 0, 3.29e1, 0, 0, 8.84, 0 };
  //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  state.conditions_[0].temperature_ = 273.15;
  state.conditions_[0].pressure_ = 1000 * 100;  // 1000 hPa
  state.conditions_[0].air_density_ = 2.7e19;
  state.custom_rate_parameters_[0][0] = 1.0e-4;
  state.custom_rate_parameters_[0][1] = 1.0e-5;
  state.custom_rate_parameters_[0][2] = 1.0e-6;
  double time_start = 0;
  double time_end = 1;

  for (auto& elem : state.variables_[0])
  {
    elem /= 10;
  }

  solver.UpdateState(state);
  auto results = solver.Solve(time_start, time_end, state);
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
TEST(RegressionChapmanODESolver, solve_10_times_smaller)
{
  micm::ChapmanODESolver solver{};
  // micm::ChapmanODESolver general{};
  TestSolve10TimesSmaller(solver);
}

void TestSolveWithRandomNumberDensities(micm::ChapmanODESolver& solver)
{
  micm::State<micm::Matrix> state = solver.GetState();
  state.conditions_[0].temperature_ = 273.15;
  state.conditions_[0].pressure_ = 1000 * 100;  // 1000 hPa
  state.conditions_[0].air_density_ = 2.7e19;
  state.custom_rate_parameters_[0][0] = 1.0e-4;
  state.custom_rate_parameters_[0][1] = 1.0e-5;
  state.custom_rate_parameters_[0][2] = 1.0e-6;
  double time_start = 0;
  double time_end = 1;

  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  std::generate(state.variables_[0].begin(), state.variables_[0].end(), [&] { return distribution(generator); });

  solver.UpdateState(state);

  auto results = solver.Solve(time_start, time_end, state);
  EXPECT_EQ(results.state_, micm::ChapmanODESolver::SolverState::Converged);
}
TEST(ChapmanMechanismHardCodedAndGeneral, solve_with_random_number_densities)
{
  micm::ChapmanODESolver solver{};
  TestSolveWithRandomNumberDensities(solver);
}
