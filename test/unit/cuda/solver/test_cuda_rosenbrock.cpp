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
void TestNormalizedErrorConst(const std::size_t number_of_grid_cells = L)
{
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  gpu_builder = GetSolver(gpu_builder);
  auto gpu_solver = gpu_builder.Build();
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
void TestNormalizedErrorDiff(const std::size_t number_of_grid_cells = L)
{
  auto gpu_builder = GpuBuilder<L>(micm::CudaRosenbrockSolverParameters::ThreeStageRosenbrockParameters());
  gpu_builder = GetSolver(gpu_builder);
  auto gpu_solver = gpu_builder.Build();
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
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  TestAlphaMinusJacobian(GpuBuilder<2>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  TestAlphaMinusJacobian(GpuBuilder<7>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  TestAlphaMinusJacobian(GpuBuilder<29>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  TestAlphaMinusJacobian(GpuBuilder<37>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  TestAlphaMinusJacobian(GpuBuilder<77>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  TestAlphaMinusJacobian(GpuBuilder<219>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  TestAlphaMinusJacobian(GpuBuilder<5599>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  TestAlphaMinusJacobian(GpuBuilder<6603>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  TestAlphaMinusJacobian(GpuBuilder<200041>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  TestAlphaMinusJacobian(GpuBuilder<421875>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  TestAlphaMinusJacobian(GpuBuilder<3395043>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);

  // number of grid cells != cuda matrix vector length
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  TestAlphaMinusJacobian(GpuBuilder<1>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);

  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 7);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 29);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 37);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 77);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 219);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 5599);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 6603);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 200041);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 421875);
  TestAlphaMinusJacobian(GpuBuilder<1109>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3395043);
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
  TestNormalizedErrorConst<1>();
  TestNormalizedErrorConst<2>();
  TestNormalizedErrorConst<7>();
  TestNormalizedErrorConst<29>();
  TestNormalizedErrorConst<37>();
  TestNormalizedErrorConst<77>();
  TestNormalizedErrorConst<219>();
  TestNormalizedErrorConst<5599>();
  TestNormalizedErrorConst<6603>();
  TestNormalizedErrorConst<200041>();
  TestNormalizedErrorConst<421875>();
  TestNormalizedErrorConst<3395043>();

  // number of grid cells != cuda matrix vector length
  TestNormalizedErrorConst<1>(2);
  TestNormalizedErrorConst<1>(7);
  TestNormalizedErrorConst<1>(29);
  TestNormalizedErrorConst<1>(37);
  TestNormalizedErrorConst<1>(77);
  TestNormalizedErrorConst<1>(219);
  TestNormalizedErrorConst<1>(5599);
  TestNormalizedErrorConst<1>(6603);
  TestNormalizedErrorConst<1>(200041);
  TestNormalizedErrorConst<1>(421875);
  TestNormalizedErrorConst<1>(3395043);

  TestNormalizedErrorConst<1109>(1);
  TestNormalizedErrorConst<1109>(2);
  TestNormalizedErrorConst<1109>(7);
  TestNormalizedErrorConst<1109>(29);
  TestNormalizedErrorConst<1109>(37);
  TestNormalizedErrorConst<1109>(77);
  TestNormalizedErrorConst<1109>(219);
  TestNormalizedErrorConst<1109>(5599);
  TestNormalizedErrorConst<1109>(6603);
  TestNormalizedErrorConst<1109>(200041);
  TestNormalizedErrorConst<1109>(421875);
  TestNormalizedErrorConst<1109>(3395043);

  /*******************************************************/
  /* tests where RMSE changes with the size of the array */
  /*******************************************************/

  // number of grid cells == cuda matrix vector length
  TestNormalizedErrorDiff<1>();
  TestNormalizedErrorDiff<2>();
  TestNormalizedErrorDiff<4>();
  TestNormalizedErrorDiff<7>();
  TestNormalizedErrorDiff<29>();
  TestNormalizedErrorDiff<37>();
  TestNormalizedErrorDiff<77>();
  TestNormalizedErrorDiff<219>();
  TestNormalizedErrorDiff<5599>();
  TestNormalizedErrorDiff<6603>();
  TestNormalizedErrorDiff<200041>();
  TestNormalizedErrorDiff<421875>();
  TestNormalizedErrorDiff<3395043>();

  // number of grid cells != cuda matrix vector length
  TestNormalizedErrorDiff<1>(2);
  TestNormalizedErrorDiff<1>(4);
  TestNormalizedErrorDiff<1>(7);
  TestNormalizedErrorDiff<1>(29);
  TestNormalizedErrorDiff<1>(37);
  TestNormalizedErrorDiff<1>(77);
  TestNormalizedErrorDiff<1>(219);
  TestNormalizedErrorDiff<1>(5599);
  TestNormalizedErrorDiff<1>(6603);
  TestNormalizedErrorDiff<1>(200041);
  TestNormalizedErrorDiff<1>(421875);
  TestNormalizedErrorDiff<1>(3395043);

  TestNormalizedErrorDiff<1109>(1);
  TestNormalizedErrorDiff<1109>(2);
  TestNormalizedErrorDiff<1109>(4);
  TestNormalizedErrorDiff<1109>(7);
  TestNormalizedErrorDiff<1109>(29);
  TestNormalizedErrorDiff<1109>(37);
  TestNormalizedErrorDiff<1109>(77);
  TestNormalizedErrorDiff<1109>(219);
  TestNormalizedErrorDiff<1109>(5599);
  TestNormalizedErrorDiff<1109>(6603);
  TestNormalizedErrorDiff<1109>(200041);
  TestNormalizedErrorDiff<1109>(421875);
  TestNormalizedErrorDiff<1109>(3395043);
}