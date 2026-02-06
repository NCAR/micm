#include "../../solver/test_rosenbrock_solver_policy.hpp"

#include <micm/GPU.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <iostream>

template<std::size_t L>
using GpuBuilder = micm::CudaSolverBuilderInPlace<micm::CudaRosenbrockSolverParameters, L>;

// In this test, all the elements in the same array are identical;
// thus the calculated RMSE should be the same no matter what the size of the array is.
template<std::size_t L>
void testNormalizedErrorConst(const std::size_t number_of_grid_cells = L)
{
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  // gpu_builder = getSolver(gpu_builder);
  // auto gpu_solver = gpu_builder.Build();
  auto gpu_solver = getSolver(std::move(gpu_builder)).Build();
  // auto gpu_solver = gpu_builder.Build();
  auto state = gpu_solver.GetState(number_of_grid_cells);
  auto& atol = state.absolute_tolerance_;
  double rtol = state.relative_tolerance_;

  auto y_old = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 1.0);
  auto y_new = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 2.0);
  auto errors = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 3.0);

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  state.SyncInputsToDevice();

  double error = gpu_solver.solver_.NormalizedError(y_old, y_new, errors, state);

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
void testNormalizedErrorDiff(const std::size_t number_of_grid_cells = L)
{
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  // gpu_builder = getSolver(gpu_builder);
  // auto gpu_solver = gpu_builder.Build();
  // gpu_builder = getSolver(std::move(gpu_builder)).Build();
  auto gpu_solver = getSolver(std::move(gpu_builder)).Build();
  auto state = gpu_solver.GetState(number_of_grid_cells);
  auto& atol = state.absolute_tolerance_;
  double rtol = state.relative_tolerance_;
  auto y_old = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 7.7);
  auto y_new = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, -13.9);
  auto errors = micm::CudaDenseMatrix<double, L>(number_of_grid_cells, state.state_size_, 81.57);

  // manually change each value of atol
  for (size_t i = 0; i < state.state_size_; ++i)
  {
    atol[i] = atol[i] * (1 + i * 0.01);
  }
  state.SetAbsoluteTolerances(atol);  // copy atol to the device

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
  state.SyncInputsToDevice();

  double computed_error = gpu_solver.solver_.NormalizedError(y_old, y_new, errors, state);

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
  // number of grid cells == cuda matrix vector length
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testAlphaMinusJacobian(GpuBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testAlphaMinusJacobian(GpuBuilder<7>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  testAlphaMinusJacobian(GpuBuilder<29>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  testAlphaMinusJacobian(GpuBuilder<37>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  testAlphaMinusJacobian(GpuBuilder<77>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  testAlphaMinusJacobian(GpuBuilder<219>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  testAlphaMinusJacobian(GpuBuilder<5599>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  testAlphaMinusJacobian(GpuBuilder<6603>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  testAlphaMinusJacobian(GpuBuilder<200041>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  testAlphaMinusJacobian(GpuBuilder<421875>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  testAlphaMinusJacobian(GpuBuilder<3395043>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);

  // number of grid cells != cuda matrix vector length
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  testAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);

  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  testAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);
}

TEST(RosenbrockSolver, CudaNormalizedError)
{
  // Based on my experience, assuming we use BLOCK_SIZE = N, it is better to test the length size (L)
  // of arrays at least within the following ranges for robustness: [L<N, N<L<2N, 2N<L<4N, 4N<L<N^2, N^2<L<N^3];
  // Here L = state_size_ * number_of_grid_cells_
  // Trying some odd and weird numbers is always helpful to reveal a potential bug.

  /***************************************************************/
  /* tests where RMSE does not change with the size of the array */
  /***************************************************************/

  // number of grid cells == cuda matrix vector length
  testNormalizedErrorConst<1>();
  testNormalizedErrorConst<2>();
  testNormalizedErrorConst<7>();
  testNormalizedErrorConst<29>();
  testNormalizedErrorConst<37>();
  testNormalizedErrorConst<77>();
  testNormalizedErrorConst<219>();
  testNormalizedErrorConst<5599>();
  testNormalizedErrorConst<6603>();
  testNormalizedErrorConst<200041>();
  testNormalizedErrorConst<421875>();
  testNormalizedErrorConst<3395043>();

  // number of grid cells != cuda matrix vector length
  testNormalizedErrorConst<1>(2);
  testNormalizedErrorConst<1>(7);
  testNormalizedErrorConst<1>(29);
  testNormalizedErrorConst<1>(37);
  testNormalizedErrorConst<1>(77);
  testNormalizedErrorConst<1>(219);
  testNormalizedErrorConst<1>(5599);
  testNormalizedErrorConst<1>(6603);
  testNormalizedErrorConst<1>(200041);
  testNormalizedErrorConst<1>(421875);
  testNormalizedErrorConst<1>(3395043);

  testNormalizedErrorConst<1109>(1);
  testNormalizedErrorConst<1109>(2);
  testNormalizedErrorConst<1109>(7);
  testNormalizedErrorConst<1109>(29);
  testNormalizedErrorConst<1109>(37);
  testNormalizedErrorConst<1109>(77);
  testNormalizedErrorConst<1109>(219);
  testNormalizedErrorConst<1109>(5599);
  testNormalizedErrorConst<1109>(6603);
  testNormalizedErrorConst<1109>(200041);
  testNormalizedErrorConst<1109>(421875);
  testNormalizedErrorConst<1109>(3395043);

  /*******************************************************/
  /* tests where RMSE changes with the size of the array */
  /*******************************************************/

  // number of grid cells == cuda matrix vector length
  testNormalizedErrorDiff<1>();
  testNormalizedErrorDiff<2>();
  testNormalizedErrorDiff<4>();
  testNormalizedErrorDiff<7>();
  testNormalizedErrorDiff<29>();
  testNormalizedErrorDiff<37>();
  testNormalizedErrorDiff<77>();
  testNormalizedErrorDiff<219>();
  testNormalizedErrorDiff<5599>();
  testNormalizedErrorDiff<6603>();
  testNormalizedErrorDiff<200041>();
  testNormalizedErrorDiff<421875>();
  testNormalizedErrorDiff<3395043>();

  // number of grid cells != cuda matrix vector length
  testNormalizedErrorDiff<1>(2);
  testNormalizedErrorDiff<1>(4);
  testNormalizedErrorDiff<1>(7);
  testNormalizedErrorDiff<1>(29);
  testNormalizedErrorDiff<1>(37);
  testNormalizedErrorDiff<1>(77);
  testNormalizedErrorDiff<1>(219);
  testNormalizedErrorDiff<1>(5599);
  testNormalizedErrorDiff<1>(6603);
  testNormalizedErrorDiff<1>(200041);
  testNormalizedErrorDiff<1>(421875);
  testNormalizedErrorDiff<1>(3395043);

  testNormalizedErrorDiff<1109>(1);
  testNormalizedErrorDiff<1109>(2);
  testNormalizedErrorDiff<1109>(4);
  testNormalizedErrorDiff<1109>(7);
  testNormalizedErrorDiff<1109>(29);
  testNormalizedErrorDiff<1109>(37);
  testNormalizedErrorDiff<1109>(77);
  testNormalizedErrorDiff<1109>(219);
  testNormalizedErrorDiff<1109>(5599);
  testNormalizedErrorDiff<1109>(6603);
  testNormalizedErrorDiff<1109>(200041);
  testNormalizedErrorDiff<1109>(421875);
  testNormalizedErrorDiff<1109>(3395043);
}