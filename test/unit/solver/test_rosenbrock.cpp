#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy> getSolver(std::size_t number_of_grid_cells)
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

  micm::Process r1 = micm::Process::create()
                         .reactants({ foo, baz })
                         .products({ yields(bar, 1), yields(quuz, 2.4) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  micm::Process r2 = micm::Process::create()
                         .reactants({ bar })
                         .products({ yields(foo, 1), yields(quz, 1.4) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 1.0e-6 }));

  micm::Process r3 = micm::Process::create().reactants({ quz }).products({}).phase(gas_phase).rate_constant(
      micm::ArrheniusRateConstant({ .A_ = 3.5e-6 }));

  return micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3 },
      micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, false));
}

template<class T>
using SparseMatrix = micm::SparseMatrix<T>;

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testAlphaMinusJacobian(std::size_t number_of_grid_cells)
{
  auto solver = getSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>(number_of_grid_cells);
  auto jacobian = solver.GetState().jacobian_;

  EXPECT_EQ(jacobian.Size(), number_of_grid_cells);
  EXPECT_EQ(jacobian[0].size(), 5);
  EXPECT_EQ(jacobian[0][0].size(), 5);
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

  solver.AlphaMinusJacobian(jacobian, 42.042);
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
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testNormalizedErrorDiff(const size_t number_of_grid_cells)
{
  auto solver = getSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>(number_of_grid_cells);
  std::vector<double> atol = solver.parameters_.absolute_tolerance_;
  double rtol = solver.parameters_.relative_tolerance_;

  auto state = solver.GetState();
  auto y_old = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, 7.7);
  auto y_new = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, -13.9);
  auto errors = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, 81.57);

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

  double computed_error = solver.NormalizedError(y_old, y_new, errors);

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

TEST(RosenbrockSolver, StandardAlphaMinusJacobian)
{
  testAlphaMinusJacobian<micm::Matrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(1);
  testAlphaMinusJacobian<micm::Matrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(2);
  testAlphaMinusJacobian<micm::Matrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(3);
  testAlphaMinusJacobian<micm::Matrix, SparseMatrix, micm::LinearSolver<double, SparseMatrix>>(4);
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;
template<class T>
using Group8VectorMatrix = micm::VectorMatrix<T, 8>;
template<class T>
using Group10VectorMatrix = micm::VectorMatrix<T, 10>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;
template<class T>
using Group8SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<8>>;
template<class T>
using Group10SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10>>;

TEST(RosenbrockSolver, DenseAlphaMinusJacobian)
{
  testAlphaMinusJacobian<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(
      1);
  testAlphaMinusJacobian<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(
      4);
  testAlphaMinusJacobian<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(
      3);
  testAlphaMinusJacobian<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(
      2);
}

TEST(RosenbrockSolver, CanSetTolerances)
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");

  foo.SetProperty("absolute tolerance", 1.0e-07);
  bar.SetProperty("absolute tolerance", 1.0e-08);

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ foo })
                         .products({ yields(bar, 1) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  for (size_t number_of_grid_cells = 1; number_of_grid_cells <= 10; ++number_of_grid_cells)
  {
    auto solver = micm::RosenbrockSolver<>(
        micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
        std::vector<micm::Process>{ r1 },
        micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells));
    EXPECT_EQ(solver.parameters_.absolute_tolerance_.size(), 2);
    EXPECT_EQ(solver.parameters_.absolute_tolerance_[0], 1.0e-07);
    EXPECT_EQ(solver.parameters_.absolute_tolerance_[1], 1.0e-08);
  }
}

TEST(RosenbrockSolver, CanOverrideTolerancesWithParameters)
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");

  foo.SetProperty("absolute tolerance", 1.0e-07);
  bar.SetProperty("absolute tolerance", 1.0e-08);

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ foo })
                         .products({ yields(bar, 1) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  for (size_t number_of_grid_cells = 1; number_of_grid_cells <= 10; ++number_of_grid_cells)
  {
    auto params = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells);
    params.absolute_tolerance_ = { 1.0e-01, 1.0e-02 };
    auto solver = micm::RosenbrockSolver<>(
        micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::vector<micm::Process>{ r1 }, params);
    EXPECT_EQ(solver.parameters_.absolute_tolerance_.size(), 2);
    EXPECT_EQ(solver.parameters_.absolute_tolerance_[0], 1.0e-01);
    EXPECT_EQ(solver.parameters_.absolute_tolerance_[1], 1.0e-02);
  }
}

TEST(RosenbrockSolver, NormalizedError) {
  // Exact fits
  testNormalizedErrorDiff<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(1);
  testNormalizedErrorDiff<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(2);
  testNormalizedErrorDiff<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(3);
  testNormalizedErrorDiff<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(4);

  // Inexact fits
  testNormalizedErrorDiff<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(1);
  testNormalizedErrorDiff<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(2);
  testNormalizedErrorDiff<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(3);
  testNormalizedErrorDiff<Group8VectorMatrix, Group8SparseVectorMatrix, micm::LinearSolver<double, Group8SparseVectorMatrix>>(5);
  testNormalizedErrorDiff<Group10VectorMatrix, Group10SparseVectorMatrix, micm::LinearSolver<double, Group10SparseVectorMatrix>>(3);
}