#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<class SolverBuilderPolicy>
SolverBuilderPolicy getSolver(SolverBuilderPolicy builder)
{
  // ---- foo  bar  baz  quz  quuz
  // foo   0    1    2    -    -
  // bar   3    4    5    -    -
  // baz   6    -    7    -    -
  // quz   -    8    -    9    -
  // quuz 10    -   11    -    12

  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");
  auto quz = micm::Species("quz");
  auto quuz = micm::Species("quuz");

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz, quz, quuz } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ foo, baz })
                         .SetProducts({ Yields(bar, 1), Yields(quuz, 2.4) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ bar })
                         .SetProducts({ Yields(foo, 1), Yields(quz, 1.4) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 1.0e-6 }));

  micm::Process r3 = micm::Process::Create().SetReactants({ quz }).SetProducts({}).SetPhase(gas_phase).SetRateConstant(
      micm::ArrheniusRateConstant({ .A_ = 3.5e-6 }));

  return builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::vector<micm::Process>{ r1, r2, r3 })
      .SetReorderState(false);
}

template<class SolverBuilderPolicy>
void testAlphaMinusJacobian(SolverBuilderPolicy builder, std::size_t number_of_grid_cells)
{
  builder = getSolver(builder);
  auto solver = builder.SetNumberOfGridCells(number_of_grid_cells).Build();
  auto jacobian = solver.GetState().jacobian_;

  EXPECT_EQ(jacobian.NumberOfBlocks(), number_of_grid_cells);
  EXPECT_EQ(jacobian.NumRows(), 5);
  EXPECT_EQ(jacobian.NumColumns(), jacobian.NumRows());
  EXPECT_EQ(jacobian[0].Size(), 5);
  EXPECT_EQ(jacobian[0][0].Size(), 5);
  EXPECT_GE(jacobian.AsVector().size(), 13 * number_of_grid_cells);

  // Initialize Jacobian matrix
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    jacobian[i_cell][0][0] = 12.2;
    jacobian[i_cell][0][1] = 24.3 * (i_cell + 2);
    jacobian[i_cell][0][2] = 42.3;
    jacobian[i_cell][1][0] = 0.43;
    jacobian[i_cell][1][1] = 23.4;
    jacobian[i_cell][1][2] = 83.4 / (i_cell + 3);
    jacobian[i_cell][2][0] = 4.74;
    jacobian[i_cell][2][2] = 6.91;
    jacobian[i_cell][3][1] = 59.1;
    jacobian[i_cell][3][3] = 83.4;
    jacobian[i_cell][4][0] = 78.5;
    jacobian[i_cell][4][2] = 53.6;
    jacobian[i_cell][4][4] = 1.0;
  }

  // Negate the Jacobian matrix
  for (auto& elem : jacobian.AsVector())
    elem = -elem;

  solver.solver_.AlphaMinusJacobian(jacobian, 42.042);
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    EXPECT_NEAR(jacobian[i_cell][0][0], 42.042 - 12.2, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][0][1], -24.3 * (i_cell + 2), 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][0][2], -42.3, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][0], -0.43, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][1], 42.042 - 23.4, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][2], -83.4 / (i_cell + 3), 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][2][0], -4.74, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][2][2], 42.042 - 6.91, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][3][1], -59.1, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][3][3], 42.042 - 83.4, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][0], -78.5, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][2], -53.6, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][4], 42.042 - 1.0, 1.0e-5);
  }
}

// In this test, the elements in the same array are different;
// thus the calculated RMSE will change when the size of the array changes.
template<class SolverBuilderPolicy>
void testNormalizedErrorDiff(SolverBuilderPolicy builder, std::size_t number_of_grid_cells)
{
  builder = getSolver(builder);
  auto solver = builder.SetNumberOfGridCells(number_of_grid_cells).Build();
  std::vector<double> atol = solver.solver_.parameters_.absolute_tolerance_;
  double rtol = solver.solver_.parameters_.relative_tolerance_;

  auto state = solver.GetState();
  using MatrixPolicy = decltype(state.variables_);
  auto y_old = MatrixPolicy(number_of_grid_cells, state.state_size_, 7.7);
  auto y_new = MatrixPolicy(number_of_grid_cells, state.state_size_, -13.9);
  auto errors = MatrixPolicy(number_of_grid_cells, state.state_size_, 81.57);

  double expected_error = 0.0;
  for (size_t i = 0; i < number_of_grid_cells; ++i)
  {
    for (size_t j = 0; j < state.state_size_; ++j)
    {
      y_old[i][j] = y_old[i][j] * i + j;
      y_new[i][j] = y_new[i][j] / (j + 1) - i;
      errors[i][j] = errors[i][j] / (i + 7) / (j + 3);
      double ymax = std::max(std::abs(y_old[i][j]), std::abs(y_new[i][j]));
      double scale = atol[j] + rtol * ymax;
      expected_error += errors[i][j] * errors[i][j] / (scale * scale);
    }
  }
  double error_min_ = 1.0e-10;
  expected_error = std::max(std::sqrt(expected_error / (number_of_grid_cells * state.state_size_)), error_min_);

  double computed_error = solver.solver_.NormalizedError(y_old, y_new, errors);

  auto relative_error =
      std::abs(computed_error - expected_error) / std::max(std::abs(computed_error), std::abs(expected_error));

  if (relative_error > 1.e-11)
  {
    std::cout << "computed_error: " << std::setprecision(12) << computed_error << std::endl;
    std::cout << "expected_error: " << std::setprecision(12) << expected_error << std::endl;
    std::cout << "relative_error: " << std::setprecision(12) << relative_error << std::endl;
    throw std::runtime_error("Fail to match computed_error and expected_error.\n");
  }
}

using StandardBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::Matrix<double>,
    micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>>;
template<std::size_t L>
using VectorBuilder = micm::CpuSolverBuilder<
    micm::RosenbrockSolverParameters,
    micm::VectorMatrix<double, L>,
    micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

TEST(RosenbrockSolver, StandardAlphaMinusJacobian)
{
  testAlphaMinusJacobian(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testAlphaMinusJacobian(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testAlphaMinusJacobian(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testAlphaMinusJacobian(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 4);
}

TEST(RosenbrockSolver, VectorAlphaMinusJacobian)
{
  testAlphaMinusJacobian(VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testAlphaMinusJacobian(VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 4);
  testAlphaMinusJacobian(VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testAlphaMinusJacobian(VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
}

TEST(RosenbrockSolver, CanSetTolerances)
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");

  foo.SetProperty("absolute tolerance", 1.0e-07);
  bar.SetProperty("absolute tolerance", 1.0e-08);

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ foo })
                         .SetProducts({ Yields(bar, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  for (size_t number_of_grid_cells = 1; number_of_grid_cells <= 10; ++number_of_grid_cells)
  {
    auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
                      .SetReactions(std::vector<micm::Process>{ r1 })
                      .SetNumberOfGridCells(number_of_grid_cells)
                      .Build();
    EXPECT_EQ(solver.solver_.parameters_.absolute_tolerance_.size(), 2);
    EXPECT_EQ(solver.solver_.parameters_.absolute_tolerance_[0], 1.0e-07);
    EXPECT_EQ(solver.solver_.parameters_.absolute_tolerance_[1], 1.0e-08);
  }
}

TEST(RosenbrockSolver, StandardNormalizedError)
{
  testNormalizedErrorDiff(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testNormalizedErrorDiff(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testNormalizedErrorDiff(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testNormalizedErrorDiff(StandardBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 4);
}

TEST(RosenbrockSolver, VectorNormalizedError)
{
  // Exact fits
  testNormalizedErrorDiff(VectorBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testNormalizedErrorDiff(VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testNormalizedErrorDiff(VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testNormalizedErrorDiff(VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 4);

  // Inexact fits
  testNormalizedErrorDiff(VectorBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testNormalizedErrorDiff(VectorBuilder<3>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testNormalizedErrorDiff(VectorBuilder<4>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testNormalizedErrorDiff(VectorBuilder<8>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5);
  testNormalizedErrorDiff(VectorBuilder<10>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
}