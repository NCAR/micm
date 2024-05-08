#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/cuda_rosenbrock.cuh>
#include <micm/solver/cuda_rosenbrock.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/cuda_dense_matrix.hpp>
#include <micm/util/cuda_sparse_matrix.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <iostream>

template<class T>
using Group1CPUVectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group20CPUVectorMatrix = micm::VectorMatrix<T, 20>;
template<class T>
using Group300CPUVectorMatrix = micm::VectorMatrix<T, 300>;
template<class T>
using Group4000CPUVectorMatrix = micm::VectorMatrix<T, 4000>;

template<class T>
using Group1GPUVectorMatrix = micm::CudaDenseMatrix<T, 1>;
template<class T>
using Group20GPUVectorMatrix = micm::CudaDenseMatrix<T, 20>;
template<class T>
using Group300GPUVectorMatrix = micm::CudaDenseMatrix<T, 300>;
template<class T>
using Group4000GPUVectorMatrix = micm::CudaDenseMatrix<T, 4000>;

template<class T>
using Group1CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group20CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<20>>;
template<class T>
using Group300CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<300>>;
template<class T>
using Group4000CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4000>>;

template<class T>
using Group1GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group20GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<20>>;
template<class T>
using Group300GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<300>>;
template<class T>
using Group4000GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<4000>>;

// the following alias works for a CudaDenserMatrix with given row and any columns
template<class T>
using Group1CudaDenseMatrix = micm::CudaDenseMatrix<T, 1>;
template<class T>
using Group2CudaDenseMatrix = micm::CudaDenseMatrix<T, 2>;
template<class T>
using Group4CudaDenseMatrix = micm::CudaDenseMatrix<T, 4>;
template<class T>
using Group7CudaDenseMatrix = micm::CudaDenseMatrix<T, 7>;
template<class T>
using Group12CudaDenseMatrix = micm::CudaDenseMatrix<T, 12>;
template<class T>
using Group16CudaDenseMatrix = micm::CudaDenseMatrix<T, 16>;
template<class T>
using Group20CudaDenseMatrix = micm::CudaDenseMatrix<T, 20>;
template<class T>
using Group5599CudaDenseMatrix = micm::CudaDenseMatrix<T, 5599>;
template<class T>
using Group6603CudaDenseMatrix = micm::CudaDenseMatrix<T, 6603>;
template<class T>
using Group200041CudaDenseMatrix = micm::CudaDenseMatrix<T, 200041>;
template<class T>
using Group421875CudaDenseMatrix = micm::CudaDenseMatrix<T, 421875>;
template<class T>
using Group3395043CudaDenseMatrix = micm::CudaDenseMatrix<T, 3395043>;

// the following alias works for a CudaSparseMatrix with given rows and any columns
template<class T>
using Group1CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group4CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;
template<class T>
using Group7CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<7>>;
template<class T>
using Group12CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<12>>;
template<class T>
using Group16CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<16>>;
template<class T>
using Group20CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<20>>;
template<class T>
using Group5599CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<5599>>;
template<class T>
using Group6603CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<6603>>;
template<class T>
using Group200041CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<200041>>;
template<class T>
using Group421875CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<421875>>;
template<class T>
using Group3395043CudaSparseMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<3395043>>;

template<class RosenbrockPolicy>
RosenbrockPolicy getSolver(std::size_t number_of_grid_cells)
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

  foo.SetProperty("absolute tolerance", 1.0e-03);
  bar.SetProperty("absolute tolerance", 1.0e-05);
  baz.SetProperty("absolute tolerance", 1.0e-07);
  quz.SetProperty("absolute tolerance", 1.0e-08);
  quuz.SetProperty("absolute tolerance", 1.0e-10);

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

  return RosenbrockPolicy(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3 },
      micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters(number_of_grid_cells, false));
}

template<
    template<class>
    class CPUMatrixPolicy,
    template<class>
    class CPUSparseMatrixPolicy,
    template<class>
    class GPUMatrixPolicy,
    template<class>
    class GPUSparseMatrixPolicy>
void testAlphaMinusJacobian(std::size_t number_of_grid_cells)
{
  auto gpu_solver = getSolver<micm::CudaRosenbrockSolver<GPUMatrixPolicy, GPUSparseMatrixPolicy>>(number_of_grid_cells);
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
  auto cpu_jacobian = gpu_jacobian;

  gpu_jacobian.CopyToDevice();
  gpu_solver.AlphaMinusJacobian(gpu_jacobian, 42.042);
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

  auto cpu_solver = getSolver<micm::RosenbrockSolver<CPUMatrixPolicy, CPUSparseMatrixPolicy>>(number_of_grid_cells);
  cpu_solver.AlphaMinusJacobian(cpu_jacobian, 42.042);

  std::vector<double> jacobian_gpu_vector = gpu_jacobian.AsVector();
  std::vector<double> jacobian_cpu_vector = cpu_jacobian.AsVector();
  for (int i = 0; i < jacobian_cpu_vector.size(); i++)
  {
    EXPECT_EQ(jacobian_cpu_vector[i], jacobian_gpu_vector[i]);
  }
}

// In this test, all the elements in the same array are identical;
// thus the calculated RMSE should be the same no matter what the size of the array is.
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testNormalizedErrorConst(const size_t number_of_grid_cells)
{
  auto gpu_solver = getSolver<micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>>(number_of_grid_cells);
  std::vector<double> atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
  auto y_old = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, 1.0);
  auto y_new = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, 2.0);
  auto errors = MatrixPolicy<double>(number_of_grid_cells, state.state_size_, 3.0);

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double error = gpu_solver.NormalizedError(y_old, y_new, errors);

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
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testNormalizedErrorDiff(const size_t number_of_grid_cells)
{
  auto gpu_solver = getSolver<micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>>(number_of_grid_cells);
  std::vector<double> atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
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

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double computed_error = gpu_solver.NormalizedError(y_old, y_new, errors);

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
  testAlphaMinusJacobian<
      Group1CPUVectorMatrix,
      Group1CPUSparseVectorMatrix,
      Group1GPUVectorMatrix,
      Group1GPUSparseVectorMatrix>(1);
  testAlphaMinusJacobian<
      Group20CPUVectorMatrix,
      Group20CPUSparseVectorMatrix,
      Group20GPUVectorMatrix,
      Group20GPUSparseVectorMatrix>(20);
  testAlphaMinusJacobian<
      Group300CPUVectorMatrix,
      Group300CPUSparseVectorMatrix,
      Group300GPUVectorMatrix,
      Group300GPUSparseVectorMatrix>(300);
  testAlphaMinusJacobian<
      Group4000CPUVectorMatrix,
      Group4000CPUSparseVectorMatrix,
      Group4000GPUVectorMatrix,
      Group4000GPUSparseVectorMatrix>(4000);
}

TEST(RosenbrockSolver, CudaNormalizedError)
{
  // Based on my experience, assuming we use BLOCK_SIZE = N, it is better to test the length size (L)
  // of arrays at least within the following ranges for robustness: [L<N, N<L<2N, 2N<L<4N, 4N<L<N^2, N^2<L<N^3];
  // Here L = state_size_ * number_of_grid_cells_
  // Trying some odd and weird numbers is always helpful to reveal a potential bug.

  // tests where RMSE does not change with the size of the array
  testNormalizedErrorConst<Group1CudaDenseMatrix, Group1CudaSparseMatrix>(1);
  testNormalizedErrorConst<Group2CudaDenseMatrix, Group2CudaSparseMatrix>(2);
  testNormalizedErrorConst<Group4CudaDenseMatrix, Group4CudaSparseMatrix>(4);
  testNormalizedErrorConst<Group7CudaDenseMatrix, Group7CudaSparseMatrix>(7);
  testNormalizedErrorConst<Group12CudaDenseMatrix, Group12CudaSparseMatrix>(12);
  testNormalizedErrorConst<Group16CudaDenseMatrix, Group16CudaSparseMatrix>(16);
  testNormalizedErrorConst<Group20CudaDenseMatrix, Group20CudaSparseMatrix>(20);
  testNormalizedErrorConst<Group5599CudaDenseMatrix, Group5599CudaSparseMatrix>(5599);
  testNormalizedErrorConst<Group6603CudaDenseMatrix, Group6603CudaSparseMatrix>(6603);
  testNormalizedErrorConst<Group200041CudaDenseMatrix, Group200041CudaSparseMatrix>(200041);
  testNormalizedErrorConst<Group421875CudaDenseMatrix, Group421875CudaSparseMatrix>(421875);
  testNormalizedErrorConst<Group3395043CudaDenseMatrix, Group3395043CudaSparseMatrix>(3395043);

  // tests where RMSE changes with the size of the array
  testNormalizedErrorDiff<Group1CudaDenseMatrix, Group1CudaSparseMatrix>(1);
  testNormalizedErrorDiff<Group2CudaDenseMatrix, Group2CudaSparseMatrix>(2);
  testNormalizedErrorDiff<Group4CudaDenseMatrix, Group4CudaSparseMatrix>(4);
  testNormalizedErrorDiff<Group7CudaDenseMatrix, Group7CudaSparseMatrix>(7);
  testNormalizedErrorDiff<Group12CudaDenseMatrix, Group12CudaSparseMatrix>(12);
  testNormalizedErrorDiff<Group16CudaDenseMatrix, Group16CudaSparseMatrix>(16);
  testNormalizedErrorDiff<Group20CudaDenseMatrix, Group20CudaSparseMatrix>(20);
  testNormalizedErrorDiff<Group5599CudaDenseMatrix, Group5599CudaSparseMatrix>(5599);
  testNormalizedErrorDiff<Group6603CudaDenseMatrix, Group6603CudaSparseMatrix>(6603);
  testNormalizedErrorDiff<Group200041CudaDenseMatrix, Group200041CudaSparseMatrix>(200041);
  testNormalizedErrorDiff<Group421875CudaDenseMatrix, Group421875CudaSparseMatrix>(421875);
  testNormalizedErrorDiff<Group3395043CudaDenseMatrix, Group3395043CudaSparseMatrix>(3395043);
}
