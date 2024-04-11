// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>

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
#include <iostream>

template<class T>
using Group1CPUVectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2CPUVectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3CPUVectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4CPUVectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1GPUVectorMatrix = micm::CudaDenseMatrix<T, 1>;
template<class T>
using Group2GPUVectorMatrix = micm::CudaDenseMatrix<T, 2>;
template<class T>
using Group3GPUVectorMatrix = micm::CudaDenseMatrix<T, 3>;
template<class T>
using Group4GPUVectorMatrix = micm::CudaDenseMatrix<T, 4>;

template<class T>
using Group1CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4CPUSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

template<class T>
using Group1GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4GPUSparseVectorMatrix = micm::CudaSparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

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

  return RosenbrockPolicy(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3 },
      micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, false));
}

template<template<class> class CPUMatrixPolicy, template<class> class CPUSparseMatrixPolicy, class CPULinearSolverPolicy,
         template<class> class GPUMatrixPolicy, template<class> class GPUSparseMatrixPolicy, class GPULinearSolverPolicy>
void testAlphaMinusJacobian(std::size_t number_of_grid_cells)
{
//  auto gpu_solver = getSolver<micm::CudaRosenbrockSolver>(number_of_grid_cells);
  auto gpu_solver = getSolver<micm::CudaRosenbrockSolver<GPUMatrixPolicy,GPUSparseMatrixPolicy>>(number_of_grid_cells);

//   auto gpu_jacobian = gpu_solver.GetState().jacobian_;
//   EXPECT_EQ(gpu_jacobian.size(), number_of_grid_cells);
//   EXPECT_EQ(gpu_jacobian[0].size(), 5);
//   EXPECT_EQ(gpu_jacobian[0][0].size(), 5);
//   auto gpu_jacobian_vec = gpu_jacobian.AsVector();
//   EXPECT_GE(gpu_jacobian_vec.size(), 13 * number_of_grid_cells);
//   gpu_jacobian_vec.assign(gpu_jacobian_vec.size(), 100.0);
//   for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
//   {
//     gpu_jacobian[i_cell][0][0] = 12.2;
//     gpu_jacobian[i_cell][0][1] = 24.3 * (i_cell + 2);
//     gpu_jacobian[i_cell][0][2] = 42.3;
//     gpu_jacobian[i_cell][1][0] = 0.43;
//     gpu_jacobian[i_cell][1][1] = 23.4;
//     gpu_jacobian[i_cell][1][2] = 83.4 / (i_cell + 3);
//     gpu_jacobian[i_cell][2][0] = 4.74;
//     gpu_jacobian[i_cell][2][2] = 6.91;
//     gpu_jacobian[i_cell][3][1] = 59.1;
//     gpu_jacobian[i_cell][3][3] = 83.4;
//     gpu_jacobian[i_cell][4][0] = 78.5;
//     gpu_jacobian[i_cell][4][2] = 53.6;
//     gpu_jacobian[i_cell][4][4] = 1.0;
//   }

//   // Negate the Jacobian matrix (-J) here
//   std::transform(gpu_jacobian_vec.cbegin(), gpu_jacobian_vec.cend(), gpu_jacobian_vec.begin(), std::negate<>{});

//   auto cpu_jacobian = gpu_jacobian;

//   gpu_jacobian.CopyToDevice();
//   gpu_solver.AlphaMinusJacobian(gpu_jacobian, 42.042);
//   gpu_jacobian.CopyToHost();

//   for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
//   {
//     EXPECT_EQ(gpu_jacobian[i_cell][0][0], 42.042 - 12.2);
//     EXPECT_EQ(gpu_jacobian[i_cell][0][1], -24.3 * (i_cell + 2));
//     EXPECT_EQ(gpu_jacobian[i_cell][0][2], -42.3);
//     EXPECT_EQ(gpu_jacobian[i_cell][1][0], -0.43);
//     EXPECT_EQ(gpu_jacobian[i_cell][1][1], 42.042 - 23.4);
//     EXPECT_EQ(gpu_jacobian[i_cell][1][2], -83.4 / (i_cell + 3));
//     EXPECT_EQ(gpu_jacobian[i_cell][2][0], -4.74);
//     EXPECT_EQ(gpu_jacobian[i_cell][2][2], 42.042 - 6.91);
//     EXPECT_EQ(gpu_jacobian[i_cell][3][1], -59.1);
//     EXPECT_EQ(gpu_jacobian[i_cell][3][3], 42.042 - 83.4);
//     EXPECT_EQ(gpu_jacobian[i_cell][4][0], -78.5);
//     EXPECT_EQ(gpu_jacobian[i_cell][4][2], -53.6);
//     EXPECT_EQ(gpu_jacobian[i_cell][4][4], 42.042 - 1.0);
//   }

//    auto cpu_solver = getSolver<
//        CPUMatrixPolicy,
//        CPUSparseMatrixPolicy,
//        micm::LinearSolver<double, CPUSparseMatrixPolicy>,
//        micm::RosenbrockSolver<CPUMatrixPolicy, CPUSparseMatrixPolicy, micm::LinearSolver<double, CPUSparseMatrixPolicy>>>(
//        number_of_grid_cells);
//   cpu_solver.AlphaMinusJacobian(cpu_jacobian, 42.042);

//   std::vector<double> jacobian_gpu_vector = gpu_jacobian.AsVector();
//   std::vector<double> jacobian_cpu_vector = cpu_jacobian.AsVector();
//   for (int i = 0; i < jacobian_cpu_vector.size(); i++)
//   {
//     EXPECT_EQ(jacobian_cpu_vector[i], jacobian_gpu_vector[i]);
//   }
}

// In this test, all the elements in the same array are identical;
// thus the calculated RMSE should be the same no matter what the size of the array is.
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testNormalizedErrorConst(const size_t num_grid_cells)
{
  auto gpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      LinearSolverPolicy,
      micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(num_grid_cells);

  double atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
  auto y_old = MatrixPolicy<double>(num_grid_cells, state.state_size_, 1.0);
  auto y_new = MatrixPolicy<double>(num_grid_cells, state.state_size_, 2.0);
  auto errors = MatrixPolicy<double>(num_grid_cells, state.state_size_, 3.0);

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double error = gpu_solver.NormalizedError(y_old, y_new, errors);

  double denom = atol + rtol * 2.0;
  // use the following function instead to avoid tiny numerical differece
  EXPECT_DOUBLE_EQ(error, std::sqrt(3.0 * 3.0 / (denom * denom)));
}

// In this test, the elements in the same array are different;
// thus the calculated RMSE will change when the size of the array changes.
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testNormalizedErrorDiff(const size_t num_grid_cells)
{
  auto gpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      LinearSolverPolicy,
      micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(num_grid_cells);

  double atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto state = gpu_solver.GetState();
  auto y_old = MatrixPolicy<double>(num_grid_cells, state.state_size_, 7.7);
  auto y_new = MatrixPolicy<double>(num_grid_cells, state.state_size_, -13.9);
  auto errors = MatrixPolicy<double>(num_grid_cells, state.state_size_, 81.57);

  double expected_error = 0.0;
  for (size_t i = 0; i < num_grid_cells; ++i)
  {
    for (size_t j = 0; j < state.state_size_; ++j)
    {
      y_old[i][j] = y_old[i][j] * i + j;
      y_new[i][j] = y_new[i][j] / (j + 1) - i;
      errors[i][j] = errors[i][j] / (i + 7) / (j + 3);
      double ymax = std::max(std::abs(y_old[i][j]), std::abs(y_new[i][j]));
      double scale = atol + rtol * ymax;
      expected_error += errors[i][j] * errors[i][j] / (scale * scale);
    }
  }
  double error_min_ = 1.0e-10;
  expected_error = std::max(std::sqrt(expected_error / (num_grid_cells * state.state_size_)), error_min_);

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
      micm::LinearSolver<double, Group1CPUSparseVectorMatrix>,
      Group1GPUVectorMatrix,
      Group1GPUSparseVectorMatrix,
      micm::CudaLinearSolver<double, Group1GPUSparseVectorMatrix>>(1);
  testAlphaMinusJacobian<
      Group2CPUVectorMatrix,
      Group2CPUSparseVectorMatrix,
      micm::LinearSolver<double, Group2CPUSparseVectorMatrix>,
      Group2GPUVectorMatrix,
      Group2GPUSparseVectorMatrix,
      micm::CudaLinearSolver<double, Group2GPUSparseVectorMatrix>>(2);
  testAlphaMinusJacobian<
      Group3CPUVectorMatrix,
      Group3CPUSparseVectorMatrix,
      micm::LinearSolver<double, Group3CPUSparseVectorMatrix>,
      Group3GPUVectorMatrix,
      Group3GPUSparseVectorMatrix,
      micm::CudaLinearSolver<double, Group3GPUSparseVectorMatrix>>(3);
  testAlphaMinusJacobian<
      Group4CPUVectorMatrix,
      Group4CPUSparseVectorMatrix,
      micm::LinearSolver<double, Group4CPUSparseVectorMatrix>,
      Group4GPUVectorMatrix,
      Group4GPUSparseVectorMatrix,
      micm::CudaLinearSolver<double, Group4GPUSparseVectorMatrix>>(4);
}

// TEST(RosenbrockSolver, CudaNormalizedError)
// {
//   // Based on my experience, assuming we use BLOCK_SIZE = N, it is better to test the length size (L)
//   // of arrays at least within the following ranges for robustness: [L<N, N<L<2N, 2N<L<4N, 4N<L<N^2, N^2<L<N^3];
//   // Here L = state_size_ * number_of_grid_cells_
//   // Trying some odd and weird numbers is always helpful to reveal a potential bug.

//   // tests where RMSE does not change with the size of the array
//   testNormalizedErrorConst<
//       Group1CudaDenseMatrix,
//       Group1CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group1CudaSparseMatrix>>(1);
//   testNormalizedErrorConst<
//       Group2CudaDenseMatrix,
//       Group2CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group2CudaSparseMatrix>>(2);
//   testNormalizedErrorConst<
//       Group4CudaDenseMatrix,
//       Group4CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group4CudaSparseMatrix>>(4);
//   testNormalizedErrorConst<
//       Group7CudaDenseMatrix,
//       Group7CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group7CudaSparseMatrix>>(7);
//   testNormalizedErrorConst<
//       Group12CudaDenseMatrix,
//       Group12CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group12CudaSparseMatrix>>(12);
//   testNormalizedErrorConst<
//       Group16CudaDenseMatrix,
//       Group16CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group16CudaSparseMatrix>>(16);
//   testNormalizedErrorConst<
//       Group20CudaDenseMatrix,
//       Group20CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group20CudaSparseMatrix>>(20);
//   testNormalizedErrorConst<
//       Group5599CudaDenseMatrix,
//       Group5599CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group5599CudaSparseMatrix>>(5599);
//   testNormalizedErrorConst<
//       Group6603CudaDenseMatrix,
//       Group6603CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group6603CudaSparseMatrix>>(6603);
//   testNormalizedErrorConst<
//       Group200041CudaDenseMatrix,
//       Group200041CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group200041CudaSparseMatrix>>(200041);
//   testNormalizedErrorConst<
//       Group421875CudaDenseMatrix,
//       Group421875CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group421875CudaSparseMatrix>>(421875);
//   testNormalizedErrorConst<
//       Group3395043CudaDenseMatrix,
//       Group3395043CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group3395043CudaSparseMatrix>>(3395043);

//   // tests where RMSE changes with the size of the array
//   testNormalizedErrorDiff<
//       Group1CudaDenseMatrix,
//       Group1CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group1CudaSparseMatrix>>(1);
//   testNormalizedErrorDiff<
//       Group2CudaDenseMatrix,
//       Group2CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group2CudaSparseMatrix>>(2);
//   testNormalizedErrorDiff<
//       Group4CudaDenseMatrix,
//       Group4CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group4CudaSparseMatrix>>(4);
//   testNormalizedErrorDiff<
//       Group7CudaDenseMatrix,
//       Group7CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group7CudaSparseMatrix>>(7);
//   testNormalizedErrorDiff<
//       Group12CudaDenseMatrix,
//       Group12CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group12CudaSparseMatrix>>(12);
//   testNormalizedErrorDiff<
//       Group16CudaDenseMatrix,
//       Group16CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group16CudaSparseMatrix>>(16);
//   testNormalizedErrorDiff<
//       Group20CudaDenseMatrix,
//       Group20CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group20CudaSparseMatrix>>(20);
//   testNormalizedErrorDiff<
//       Group5599CudaDenseMatrix,
//       Group5599CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group5599CudaSparseMatrix>>(5599);
//   testNormalizedErrorDiff<
//       Group6603CudaDenseMatrix,
//       Group6603CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group6603CudaSparseMatrix>>(6603);
//   testNormalizedErrorDiff<
//       Group200041CudaDenseMatrix,
//       Group200041CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group200041CudaSparseMatrix>>(200041);
//   testNormalizedErrorDiff<
//       Group421875CudaDenseMatrix,
//       Group421875CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group421875CudaSparseMatrix>>(421875);
//   testNormalizedErrorDiff<
//       Group3395043CudaDenseMatrix,
//       Group3395043CudaSparseMatrix,
//       micm::CudaLinearSolver<double, Group3395043CudaSparseMatrix>>(3395043);
// }
