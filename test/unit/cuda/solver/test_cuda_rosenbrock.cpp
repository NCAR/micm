#include "../../solver/test_rosenbrock_solver_policy.hpp"

#include <micm/cuda/solver/cuda_rosenbrock.cuh>
#include <micm/cuda/solver/cuda_rosenbrock.hpp>
#include <micm/cuda/solver/cuda_solver_builder.hpp>
#include <micm/cuda/solver/cuda_solver_parameters.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <iostream>

template<std::size_t L>
using GpuBuilder = micm::CudaSolverBuilder<micm::CudaRosenbrockSolverParameters, L>;

template<std::size_t L>
void testAlphaMinusJacobian()
{
  std::size_t number_of_grid_cells = L;
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  gpu_builder = getSolver(gpu_builder);
  auto gpu_solver = gpu_builder.SetNumberOfGridCells(number_of_grid_cells).Build();
  auto cpu_builder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  cpu_builder = getSolver(cpu_builder);
  auto cpu_solver = cpu_builder.SetNumberOfGridCells(number_of_grid_cells).Build();

  auto gpu_jacobian = gpu_solver.GetState().jacobian_;
  EXPECT_EQ(gpu_jacobian.NumberOfBlocks(), number_of_grid_cells);
  EXPECT_EQ(gpu_jacobian.NumRows(), 5);
  EXPECT_EQ(gpu_jacobian.NumColumns(), gpu_jacobian.NumRows());
  EXPECT_EQ(gpu_jacobian[0].Size(), 5);
  EXPECT_EQ(gpu_jacobian[0][0].Size(), 5);
  auto& gpu_jacobian_vec = gpu_jacobian.AsVector();
  EXPECT_GE(gpu_jacobian_vec.size(), 13 * number_of_grid_cells);

  gpu_jacobian_vec.assign(gpu_jacobian_vec.size(), 100.0);
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    gpu_jacobian[i_cell][0][0] = 12.2;
    gpu_jacobian[i_cell][0][1] = 24.3 * (i_cell + 2);
    gpu_jacobian[i_cell][0][2] = 42.3;
    gpu_jacobian[i_cell][1][0] = 0.43;
    gpu_jacobian[i_cell][1][1] = 23.4;
    gpu_jacobian[i_cell][1][2] = 83.4 / (i_cell + 3);
    gpu_jacobian[i_cell][2][0] = 4.74;
    gpu_jacobian[i_cell][2][2] = 6.91;
    gpu_jacobian[i_cell][3][1] = 59.1;
    gpu_jacobian[i_cell][3][3] = 83.4;
    gpu_jacobian[i_cell][4][0] = 78.5;
    gpu_jacobian[i_cell][4][2] = 53.6;
    gpu_jacobian[i_cell][4][4] = 1.0;
  }

  // Negate the Jacobian matrix (-J) here
  std::transform(gpu_jacobian_vec.cbegin(), gpu_jacobian_vec.cend(), gpu_jacobian_vec.begin(), std::negate<>{});

  auto cpu_jacobian = cpu_solver.GetState().jacobian_;
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i = 0; i < 5; ++i)
    {
      for (std::size_t j = 0; j < 5; ++j)
      {
        if (!cpu_jacobian.IsZero(i, j))
          cpu_jacobian[i_cell][i][j] = gpu_jacobian[i_cell][i][j];
      }
    }
  }

  gpu_jacobian.CopyToDevice();
  gpu_solver.solver_.AlphaMinusJacobian(gpu_jacobian, 42.042);
  gpu_jacobian.CopyToHost();

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    EXPECT_EQ(gpu_jacobian[i_cell][0][0], 42.042 - 12.2);
    EXPECT_EQ(gpu_jacobian[i_cell][0][1], -24.3 * (i_cell + 2));
    EXPECT_EQ(gpu_jacobian[i_cell][0][2], -42.3);
    EXPECT_EQ(gpu_jacobian[i_cell][1][0], -0.43);
    EXPECT_EQ(gpu_jacobian[i_cell][1][1], 42.042 - 23.4);
    EXPECT_EQ(gpu_jacobian[i_cell][1][2], -83.4 / (i_cell + 3));
    EXPECT_EQ(gpu_jacobian[i_cell][2][0], -4.74);
    EXPECT_EQ(gpu_jacobian[i_cell][2][2], 42.042 - 6.91);
    EXPECT_EQ(gpu_jacobian[i_cell][3][1], -59.1);
    EXPECT_EQ(gpu_jacobian[i_cell][3][3], 42.042 - 83.4);
    EXPECT_EQ(gpu_jacobian[i_cell][4][0], -78.5);
    EXPECT_EQ(gpu_jacobian[i_cell][4][2], -53.6);
    EXPECT_EQ(gpu_jacobian[i_cell][4][4], 42.042 - 1.0);
  }

  cpu_solver.solver_.AlphaMinusJacobian(cpu_jacobian, 42.042);

  // Compare the results
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    for (std::size_t i = 0; i < 5; ++i)
    {
      for (std::size_t j = 0; j < 5; ++j)
      {
        if (!cpu_jacobian.IsZero(i, j))
          EXPECT_EQ(cpu_jacobian[i_cell][i][j], gpu_jacobian[i_cell][i][j]);
      }
    }
  }
}

// In this test, all the elements in the same array are identical;
// thus the calculated RMSE should be the same no matter what the size of the array is.
template<std::size_t L>
void testNormalizedErrorConst()
{
  std::size_t number_of_grid_cells = L;
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  gpu_builder = getSolver(gpu_builder);
  auto gpu_solver = gpu_builder.SetNumberOfGridCells(number_of_grid_cells).Build();

  std::vector<double> atol = gpu_solver.solver_.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.solver_.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
  auto y_old = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 1.0);
  auto y_new = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 2.0);
  auto errors = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 3.0);

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double error = gpu_solver.solver_.NormalizedError(y_old, y_new, errors);

  auto expected_error = 0.0;
  for (size_t i = 0; i < state.state_size_; ++i)
  {
    double ymax = std::max(std::abs(y_old[0][i]), std::abs(y_new[0][i]));
    double scale = atol[i] + rtol * ymax;
    expected_error += errors[0][i] * errors[0][i] / (scale * scale);
  }
  double error_min_ = 1.0e-10;
  expected_error = std::max(std::sqrt(expected_error / state.state_size_), error_min_);

  // use the following function instead to avoid tiny numerical differece
  auto relative_error = std::abs(error - expected_error) / std::max(std::abs(error), std::abs(expected_error));
  if (relative_error > 1.e-14)
  {
    std::cout << "error: " << std::setprecision(12) << error << std::endl;
    std::cout << "expected_error: " << std::setprecision(12) << expected_error << std::endl;
    std::cout << "relative_error: " << std::setprecision(12) << relative_error << std::endl;
    throw std::runtime_error("Fail to match error and expected_error.\n");
  }
}

// In this test, the elements in the same array are different;
// thus the calculated RMSE will change when the size of the array changes.
template<std::size_t L>
void testNormalizedErrorDiff()
{
  std::size_t number_of_grid_cells = L;
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  gpu_builder = getSolver(gpu_builder);
  auto gpu_solver = gpu_builder.SetNumberOfGridCells(number_of_grid_cells).Build();

  std::vector<double> atol = gpu_solver.solver_.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.solver_.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
  auto y_old = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 7.7);
  auto y_new = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, -13.9);
  auto errors = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 81.57);

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

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double computed_error = gpu_solver.solver_.NormalizedError(y_old, y_new, errors);

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

TEST(RosenbrockSolver, DenseAlphaMinusJacobian)
{
  testAlphaMinusJacobian<1>();
  testAlphaMinusJacobian<20>();
  testAlphaMinusJacobian<300>();
  testAlphaMinusJacobian<4000>();
}

TEST(RosenbrockSolver, CudaNormalizedError)
{
  // Based on my experience, assuming we use BLOCK_SIZE = N, it is better to test the length size (L)
  // of arrays at least within the following ranges for robustness: [L<N, N<L<2N, 2N<L<4N, 4N<L<N^2, N^2<L<N^3];
  // Here L = state_size_ * number_of_grid_cells_
  // Trying some odd and weird numbers is always helpful to reveal a potential bug.

  // tests where RMSE does not change with the size of the array
  testNormalizedErrorConst<1>();
  testNormalizedErrorConst<2>();
  testNormalizedErrorConst<4>();
  testNormalizedErrorConst<7>();
  testNormalizedErrorConst<12>();
  testNormalizedErrorConst<16>();
  testNormalizedErrorConst<20>();
  testNormalizedErrorConst<5599>();
  testNormalizedErrorConst<6603>();
  testNormalizedErrorConst<200041>();
  testNormalizedErrorConst<421875>();
  testNormalizedErrorConst<3395043>();

  // tests where RMSE changes with the size of the array
  testNormalizedErrorDiff<1>();
  testNormalizedErrorDiff<2>();
  testNormalizedErrorDiff<4>();
  testNormalizedErrorDiff<7>();
  testNormalizedErrorDiff<12>();
  testNormalizedErrorDiff<16>();
  testNormalizedErrorDiff<20>();
  testNormalizedErrorDiff<5599>();
  testNormalizedErrorDiff<6603>();
  testNormalizedErrorDiff<200041>();
  testNormalizedErrorDiff<421875>();
  testNormalizedErrorDiff<3395043>();
}

TEST(RosenbrockSolver, SingularSystemZeroInBottomRightOfU)
{
  auto params = micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.check_singularity_ = true;
  auto vector = GpuBuilder<4>(params);

  auto vector_solver = getSingularSystemZeroInBottomRightOfU(vector).SetNumberOfGridCells(4).Build();

  auto vector_state = vector_solver.GetState();

  double k1 = -2;
  double k2 = 1.0;

  vector_state.SetCustomRateParameter("r1", { k1, k1, k1, k1 });
  vector_state.SetCustomRateParameter("r2", { k2, k2, k2, k2 });

  vector_state.variables_[0] = { 1.0, 1.0 };
  vector_state.variables_[1] = { 1.0, 1.0 };
  vector_state.variables_[2] = { 1.0, 1.0 };
  vector_state.variables_[3] = { 1.0, 1.0 };

  // to get a jacobian with an LU factorization that contains a zero on the diagonal
  // of U, we need det(alpha * I - jacobian) = 0
  // for the system above, that means we have to have alpha + k1 + k2 = 0
  // in this case, one of the reaction rates will be negative but it's good enough to
  // test the singularity check
  // alpha is 1 / (H * gamma), where H is the time step and gamma is the gamma value from
  // the rosenbrock paramters
  // so H needs to be 1 / ( (-k1 - k2) * gamma)
  // since H is positive we need -k1 -k2 to be positive, hence the smaller, negative value for k1
  double H = 1 / ((-k1 - k2) * params.gamma_[0]);
  vector_solver.solver_.parameters_.h_start_ = H;

  vector_solver.CalculateRateConstants(vector_state);
  vector_state.SyncInputsToDevice();

  auto vector_result = vector_solver.Solve(2 * H, vector_state);
  vector_state.SyncOutputsToHost();
  EXPECT_NE(vector_result.stats_.singular_, 0);
}

TEST(RosenbrockSolver, SingularSystemZeroAlongDiagonalNotBottomRight)
{
  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  double k1 = -1.0;
  double k2 = -1.0;
  double k3 = 1.0;

  // to get a jacobian with an LU factorization that contains a zero on the diagonal
  // of U, we need det(alpha * I - jacobian) = 0
  // for the system above, that means we have to set alpha = -k1, or alpha=-k2, or alpha=k3
  double H = 1 / (-k1 * params.gamma_[0]);

  params.check_singularity_ = true;
  params.h_start_ = H;

  auto vector = GpuBuilder<4>(params);

  auto vector_solver = getSolverForSingularSystemOnDiagonal(vector).SetNumberOfGridCells(4).Build();

  auto vector_state = vector_solver.GetState();

  vector_state.SetCustomRateParameter("r1", { k1, k1, k1, k1 });
  vector_state.SetCustomRateParameter("r2", { k2, k2, k2, k2 });
  vector_state.SetCustomRateParameter("r3", { k3, k3, k3, k3 });

  vector_state.variables_[0] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[1] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[2] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[3] = { 1.0, 1.0, 1.0 };

  vector_solver.CalculateRateConstants(vector_state);
  vector_state.SyncInputsToDevice();

  auto vector_result = vector_solver.Solve(2 * H, vector_state);
  vector_state.SyncOutputsToHost();
  EXPECT_NE(vector_result.stats_.singular_, 0);
}