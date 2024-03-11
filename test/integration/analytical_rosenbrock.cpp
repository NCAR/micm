#include <gtest/gtest.h>

#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>

#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"
#include "e5.hpp"
#include "hires.hpp"
#include "oregonator.hpp"

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

TEST(AnalyticalExamples, Troe)
{
  test_analytical_troe<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Photolysis)
{
  test_analytical_photolysis<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Tunneling)
{
  test_analytical_tunneling<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Arrhenius)
{
  test_analytical_arrhenius<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Branched)
{
  test_analytical_branched<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Robertson)
{
  test_analytical_robertson<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, SurfaceRxn)
{
  test_analytical_surface_rxn<micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>>(
      [](const micm::System& s,
         const std::vector<micm::Process>& p) -> micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>
      {
        return micm::RosenbrockSolver<micm::Matrix, SparseMatrixTest>{
          s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamples, Oregonator)
{
  /*
   * I think these are the equations, but I'm really not sure. I don't know how this translates to the jacobian
   * and forcing functions used by the ODE book: https://www.unige.ch/~hairer/testset/stiff/orego/equation.f
   * A+Y -> X+P
   * X+Y -> 2P
   * A+X -> 2X+2Z
   * 2X -> A+P
   * B+Z -> 1/2fY
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ b })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .phase(gas_phase);

  auto params = micm::RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters();
  params.relative_tolerance_ = 1e-4;
  params.absolute_tolerance_ = 1e-6 * params.relative_tolerance_;
  Oregonator<micm::Matrix, SparseMatrixTest> solver(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::vector<micm::Process>{ r1, r2, r3 }, params);

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

  auto state = solver.GetState();

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
      state.variables_[0] = result.result_.AsVector();
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

TEST(AnalyticalExamples, Oregonator2)
{
  /* Equations derived from the forcing function here: https://www.unige.ch/~hairer/testset/stiff/orego/equation.f
   * a + b -> ( 1 - (1/77.27)^2 ) b    k = 77.27
   * c -> ( 1 / (0.161 * 77.27) ) b    k = 0.161
   * b -> ( 77.27 )^2 a                k = 1/77.27
   * a -> 2 a + ( 0.161/77.27 ) c      k = 77.27
   * a + a -> NULL                     k = 77.27 * 8.375e-6
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ a, b })
                         .products({ yields(b, 1 - std::pow((1 / 77.27), 2)) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ c })
                         .products({ yields(b, 1 / (0.161 * 77.27)) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create()
                         .reactants({ b })
                         .products({ yields(a, std::pow(77.27, 2)) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .phase(gas_phase);

  micm::Process r4 = micm::Process::create()
                         .reactants({ a })
                         .products({ yields(a, 2), yields(c, 0.161 / 77.27) })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .phase(gas_phase);

  micm::Process r5 = micm::Process::create()
                         .reactants({ a, a })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r5" }))
                         .phase(gas_phase);

  auto params = micm::RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters();
  params.relative_tolerance_ = 1e-4;
  params.absolute_tolerance_ = 1e-6 * params.relative_tolerance_;
  Oregonator<micm::Matrix, SparseMatrixTest> solver(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3, r4, r5 },
      params);

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

  auto state = solver.GetState();

  double k1 = 77.27;
  double k2 = 0.161;
  double k3 = 1 / 77.27;
  double k4 = 77.27;
  double k5 = 77.27 * 8.375e-6;

  state.SetCustomRateParameter("r1", k1);
  state.SetCustomRateParameter("r2", k2);
  state.SetCustomRateParameter("r3", k3);
  state.SetCustomRateParameter("r4", k4);
  state.SetCustomRateParameter("r5", k5);

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
      state.variables_[0] = result.result_.AsVector();
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
   * No idea what these equations are
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto y1 = micm::Species("y1");
  auto y2 = micm::Species("y2");
  auto y3 = micm::Species("y3");
  auto y4 = micm::Species("y4");
  auto y5 = micm::Species("y5");
  auto y6 = micm::Species("y6");
  auto y7 = micm::Species("y7");
  auto y8 = micm::Species("y8");

  micm::Phase gas_phase{ std::vector<micm::Species>{ y1, y2, y3, y4, y5, y6, y7, y8 } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ y1 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .phase(gas_phase);
  micm::Process r2 = micm::Process::create()
                         .reactants({ y2 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .phase(gas_phase);
  micm::Process r3 = micm::Process::create()
                         .reactants({ y3 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .phase(gas_phase);
  micm::Process r4 = micm::Process::create()
                         .reactants({ y4 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .phase(gas_phase);
  micm::Process r5 = micm::Process::create()
                         .reactants({ y5 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r5" }))
                         .phase(gas_phase);
  micm::Process r6 = micm::Process::create()
                         .reactants({ y6 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r6" }))
                         .phase(gas_phase);
  micm::Process r7 = micm::Process::create()
                         .reactants({ y7 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r7" }))
                         .phase(gas_phase);
  micm::Process r8 = micm::Process::create()
                         .reactants({ y8 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r8" }))
                         .phase(gas_phase);

  auto params = micm::RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters();
  params.relative_tolerance_ = 1e-3;
  params.absolute_tolerance_ = params.relative_tolerance_ * 1e-4;
  HIRES<micm::Matrix, SparseMatrixTest> solver(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7, r8 },
      params);

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

  auto state = solver.GetState();

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
      state.variables_[0] = result.result_.AsVector();
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
   * No idea what these equations are
   *
   * this problem is described in
   * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
   * edition. ed. Springer, Berlin ; New York. Page 3
   *
   * solutions are provided here
   * https://www.unige.ch/~hairer/testset/testset.html
   */

  auto y1 = micm::Species("y1");
  auto y2 = micm::Species("y2");
  auto y3 = micm::Species("y3");
  auto y4 = micm::Species("y4");

  micm::Phase gas_phase{ std::vector<micm::Species>{ y1, y2, y3, y4 } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ y1 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .phase(gas_phase);
  micm::Process r2 = micm::Process::create()
                         .reactants({ y2 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .phase(gas_phase);
  micm::Process r3 = micm::Process::create()
                         .reactants({ y3 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .phase(gas_phase);
  micm::Process r4 = micm::Process::create()
                         .reactants({ y4 })
                         .rate_constant(micm::UserDefinedRateConstant({ .label_ = "r4" }))
                         .phase(gas_phase);

  auto params = micm::RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters();
  params.relative_tolerance_ = 1e-2;
  params.absolute_tolerance_ = 1.7e-24;
  E5<micm::Matrix, SparseMatrixTest> solver(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::vector<micm::Process>{ r1, r2, r3, r4 }, params);

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

  auto state = solver.GetState();

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
      state.variables_[0] = result.result_.AsVector();
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
