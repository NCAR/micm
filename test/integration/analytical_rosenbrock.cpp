#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include "oregonator.hpp"
#include "e5.hpp"
#include "hires.hpp"

#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<double>;

TEST(AnalyticalExamples, Troe)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_troe(builder);
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_troe(builder);
}

TEST(AnalyticalExamples, Photolysis)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_photolysis(builder);
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_photolysis(builder);
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_ternary_chemical_activation(builder);
}

TEST(AnalyticalExamples, Tunneling)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_tunneling(builder);
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_tunneling(builder);
}

TEST(AnalyticalExamples, Arrhenius)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_arrhenius(builder);
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_arrhenius(builder);
}

TEST(AnalyticalExamples, Branched)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_branched(builder);
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_stiff_branched(builder);
}

TEST(AnalyticalExamples, Robertson)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_robertson(builder);
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  auto builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  test_analytical_surface_rxn(builder);
}

using LinearSolverTest = micm::LinearSolver<SparseMatrixTest, micm::LuDecomposition>;
template<class RatesPolicy>
using RosenbrockTest = micm::RosenbrockSolver<RatesPolicy, LinearSolverTest>;

TEST(AnalyticalExamples, Oregonator)
{
  /*
   * This problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using OregonatorTest = Oregonator<micm::Matrix<double>, SparseMatrixTest>;
  auto solver = OregonatorTest::template CreateSolver<RosenbrockTest<OregonatorTest>, LinearSolverTest>(params, 1);
  
  double end = 360;
  double time_step = 30;
  size_t N = static_cast<size_t>(end / time_step);

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(3));
  std::vector<std::vector<double>> analytical_concentrations(13, std::vector<double>(3));

  model_concentrations[0] = { 1, 2, 3 };

  analytical_concentrations = {
    { 1, 2, 3 },
    { 0.1000661467180497E+01, 0.1512778937348249E+04, 0.1035854312767229E+05 },
    { 0.1000874625199626E+01, 0.1144336972384497E+04, 0.8372149966624639E+02 },
    { 0.1001890368438751E+01, 0.5299926232295553E+03, 0.1662279579042420E+01 },
    { 0.1004118022612645E+01, 0.2438326079910346E+03, 0.1008822224048647E+01 },
    { 0.1008995416634061E+01, 0.1121664388662539E+03, 0.1007783229065319E+01 },
    { 0.1019763472537298E+01, 0.5159761322947535E+02, 0.1016985778956374E+01 },
    { 0.1043985088527474E+01, 0.2373442027531524E+02, 0.1037691843544522E+01 },
    { 0.1100849071667922E+01, 0.1091533805469020E+02, 0.1085831969810860E+01 },
    { 0.1249102130020572E+01, 0.5013945178605446E+01, 0.1208326626237875E+01 },
    { 0.1779724751937019E+01, 0.2281852385542403E+01, 0.1613754023671725E+01 },
    { 0.1000889326903503E+01, 0.1125438585746596E+04, 0.1641049483777168E+05 },
    { 0.1000814870318523E+01, 0.1228178521549889E+04, 0.1320554942846513E+03 },
  };

  auto state = solver.rates_.GetState();

  state.variables_[0] = model_concentrations[0];

  std::vector<double> times;
  times.push_back(0);
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      state.variables_[0] = state.variables_.AsVector();
      actual_solve += result.final_time_;
    }
    model_concentrations[i_time + 1] = state.variables_[0];
  }

  std::vector<std::string> header = { "time", "A", "B", "C" };
  writeCSV("model_concentrations.csv", header, model_concentrations, times);
  std::vector<double> an_times;
  an_times.push_back(0);
  for (int i = 1; i <= 12; ++i)
  {
    an_times.push_back(30 * i);
  }
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations, an_times);

  auto map = state.variable_map_;

  size_t _a = map.at("A");
  size_t _b = map.at("B");
  size_t _c = map.at("C");

  double tol = 1e-3;
  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    double rel_diff = relative_difference(model_concentrations[i][_a], analytical_concentrations[i][0]);
    EXPECT_TRUE(rel_diff < tol) << "Arrays differ at index (" << i << ", " << 0 << ")";
    rel_diff = relative_difference(model_concentrations[i][_b], analytical_concentrations[i][1]);
    EXPECT_TRUE(rel_diff < tol) << "Arrays differ at index (" << i << ", " << 1 << ")";
    rel_diff = relative_difference(model_concentrations[i][_c], analytical_concentrations[i][2]);
    EXPECT_TRUE(rel_diff < tol) << "Arrays differ at index (" << i << ", " << 2 << ")";
  }
}

TEST(AnalyticalExamples, HIRES)
{
  /*
   * This problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using HIRESTest = HIRES<micm::Matrix<double>, SparseMatrixTest>;
  auto solver = HIRESTest::CreateSolver<RosenbrockTest<HIRESTest>, LinearSolverTest>(params, 1);
  
  size_t N = 2;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(8));
  std::vector<std::vector<double>> analytical_concentrations(3, std::vector<double>(8));

  model_concentrations[0] = { 1, 0, 0, 0, 0, 0, 0, 0.0057 };

  analytical_concentrations = {
    { 1, 0, 0, 0, 0, 0, 0, 0.0057 },
    { 0.000737131257332567,
      0.000144248572631618,
      0.000058887297409676,
      0.001175651343283149,
      0.002386356198831330,
      0.006238968252742796,
      0.002849998395185769,
      0.002850001604814231 },
    { 0.000670305503581864,
      0.000130996846986347,
      0.000046862231597733,
      0.001044668020551705,
      0.000594883830951485,
      0.001399628833942774,
      0.001014492757718480,
      0.004685507242281520 },
  };

  auto state = solver.rates_.GetState();

  state.variables_[0] = model_concentrations[0];

  std::vector<double> times;
  times.push_back(0);
  double time_step = 321.8122;
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      state.variables_[0] = state.variables_.AsVector();
      actual_solve += result.final_time_;
    }
    model_concentrations[i_time + 1] = state.variables_[0];
    time_step += 100;
  }

  std::vector<std::string> header = { "time", "y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8" };
  writeCSV("model_concentrations.csv", header, model_concentrations, times);
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations, times);

  double tol = 1e-5;
  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    for (size_t j = 0; j < model_concentrations[0].size(); ++j)
    {
      double rel_diff = relative_difference(model_concentrations[i][j], analytical_concentrations[i][j]);
      EXPECT_NEAR(model_concentrations[i][j], analytical_concentrations[i][j], tol);
    }
  }
}

TEST(AnalyticalExamples, E5)
{
  /*
   * This problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto params = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters();
  using E5Test = E5<micm::Matrix<double>, SparseMatrixTest>;
  auto solver = E5Test::CreateSolver<RosenbrockTest<E5Test>, LinearSolverTest>(params, 1);

  size_t N = 7;

  std::vector<std::vector<double>> model_concentrations(N + 1, std::vector<double>(4));
  std::vector<std::vector<double>> analytical_concentrations(N + 1, std::vector<double>(4));

  model_concentrations[0] = { 1.76e-3, 0, 0, 0 };

  analytical_concentrations = {
    { 1.76e-3, 0, 0, 0 },
    { 1.7599259497677897058e-003, 1.3846281519376516449e-011, 7.6370038530073911180e-013, 1.3082581134075777338e-011 },
    { 1.6180769999072942552e-003, 1.3822370304983735443e-010, 8.2515735006838336088e-012, 1.2997212954915352082e-010 },
    { 7.4813208224292220114e-006, 2.3734781561205975019e-012, 2.2123586689581663654e-012, 1.6111948716243113653e-013 },
    { 4.7150333630401632232e-010, 1.8188895860807021729e-014, 1.8188812376786725407e-014, 8.3484020296321693074e-020 },
    { 3.1317148329356996037e-014, 1.4840957952870064294e-016, 1.4840957948345691466e-016, 4.5243728279782625194e-026 },
    { 3.8139035189787091771e-049, 1.0192582567660293322e-020, 1.0192582567660293322e-020, 3.7844935507486221171e-065 },
    { 0.0000000000000000000e-000, 8.8612334976263783420e-023, 8.8612334976263783421e-023, 0.0000000000000000000e-000 }
  };

  auto state = solver.rates_.GetState();

  state.variables_[0] = model_concentrations[0];

  std::vector<double> times;
  times.push_back(0);
  double time_step = 10;
  for (size_t i_time = 0; i_time < N; ++i_time)
  {
    double solve_time = time_step + i_time * time_step;
    times.push_back(solve_time);
    // Model results
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      auto result = solver.Solve(time_step - actual_solve, state);
      state.variables_[0] = state.variables_.AsVector();
      actual_solve += result.final_time_;
    }
    model_concentrations[i_time + 1] = state.variables_[0];
    time_step *= 100;
  }

  std::vector<std::string> header = { "time", "y1", "y2", "y3", "y4" };
  writeCSV("model_concentrations.csv", header, model_concentrations, times);
  writeCSV("analytical_concentrations.csv", header, analytical_concentrations, times);

  double tol = 1e-5;
  for (size_t i = 0; i < model_concentrations.size(); ++i)
  {
    for (size_t j = 0; j < model_concentrations[0].size(); ++j)
    {
      double rel_diff = relative_difference(model_concentrations[i][j], analytical_concentrations[i][j]);
      EXPECT_NEAR(model_concentrations[i][j], analytical_concentrations[i][j], tol)
          << "difference at (" << i << ", " << j << ")";
    }
  }
}