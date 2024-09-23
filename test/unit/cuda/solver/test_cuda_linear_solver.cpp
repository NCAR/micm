#include "../../solver/test_linear_solver_policy.hpp"

#include <micm/cuda/solver/cuda_linear_solver.hpp>
#include <micm/cuda/solver/cuda_lu_decomposition.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <random>

using FloatingPointType = double;

using Group10000VectorMatrix = micm::VectorMatrix<FloatingPointType, 10000>;

using Group1CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1>;
using Group20CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 20>;
using Group300CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 300>;
using Group4000CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 4000>;
using Group10000CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 10000>;

using Group10000SparseVectorMatrix = micm::SparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<10000>>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group20CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<20>>;
using Group300CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<300>>;
using Group4000CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<4000>>;
using Group10000CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<10000>>;

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
std::vector<double> linearSolverGenerator(std::size_t number_of_blocks)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, 10, 0.0);
  MatrixPolicy x(number_of_blocks, 10, 100.0);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(A);
  CopyToDeviceDense<MatrixPolicy>(b);
  CopyToDeviceDense<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 1.0e-30);
  std::pair<SparseMatrixPolicy, SparseMatrixPolicy> lu =
      micm::LuDecomposition::GetLUMatrices<SparseMatrixPolicy>(A, 1.0e-30);
  SparseMatrixPolicy lower_matrix = std::move(lu.first);
  SparseMatrixPolicy upper_matrix = std::move(lu.second);

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(lower_matrix);
  CopyToDeviceSparse<SparseMatrixPolicy>(upper_matrix);

  bool singular{ false };
  solver.Factor(A, lower_matrix, upper_matrix, singular);
  solver.template Solve<MatrixPolicy>(x, lower_matrix, upper_matrix);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<MatrixPolicy>(x);

  return x.AsVector();
}

// bit to bit variation between CPU and GPU result with randomMatrixVectorOrdering
void verify_gpu_against_cpu()
{
  std::vector<double> cpu_x = linearSolverGenerator<
      Group10000VectorMatrix,
      Group10000SparseVectorMatrix,
      micm::LinearSolver<Group10000SparseVectorMatrix>>(10000);

  std::vector<double> gpu_x = linearSolverGenerator<
      Group10000CudaDenseMatrix,
      Group10000CudaSparseMatrix,
      micm::CudaLinearSolver<Group10000CudaSparseMatrix>>(10000);

  for (int i = 0; i < cpu_x.size(); i++)
  {
    auto cpu_ele = cpu_x[i];
    auto gpu_ele = gpu_x[i];
    EXPECT_LT(std::abs((cpu_ele - gpu_ele) / cpu_ele), 1.0e-6);
  }
}

// TEST(CudaLinearSolver, DenseMatrixVectorOrderingPolicy)
// {
//   testDenseMatrix<
//       Group1CudaDenseMatrix,
//       Group1CudaSparseMatrix,
//       micm::CudaLinearSolver<Group1CudaSparseMatrix, micm::CudaLuDecomposition>>();
// }

// TEST(CudaLinearSolver, RandomMatrixVectorOrderingPolicy)
// {
//   testRandomMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<Group1CudaSparseMatrix>>(1);
//   testRandomMatrix<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolver<Group20CudaSparseMatrix>>(20);
//   testRandomMatrix<Group300CudaDenseMatrix, Group300CudaSparseMatrix, micm::CudaLinearSolver<Group300CudaSparseMatrix>>(300);
//   testRandomMatrix<Group4000CudaDenseMatrix, Group4000CudaSparseMatrix, micm::CudaLinearSolver<Group4000CudaSparseMatrix>>(
//       4000);
// }

// TEST(CudaLinearSolver, DiagonalMatrixVectorOrderingPolicy)
// {
//   testDiagonalMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<Group1CudaSparseMatrix>>(1);
//   testDiagonalMatrix<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolver<Group20CudaSparseMatrix>>(20);
//   testDiagonalMatrix<Group300CudaDenseMatrix, Group300CudaSparseMatrix, micm::CudaLinearSolver<Group300CudaSparseMatrix>>(
//       300);
//   testDiagonalMatrix<Group4000CudaDenseMatrix, Group4000CudaSparseMatrix, micm::CudaLinearSolver<Group4000CudaSparseMatrix>>(
//       4000);
// }

// TEST(CudaLinearSolver, RandomMatrixVectorOrderingForGPU)
// {
//   verify_gpu_against_cpu();
// }

TEST(CudaLinearSolver, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  testExtremeInitialValue<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<Group1CudaSparseMatrix>>(1, INFINITY);
  // for(auto initial_value : initial_values)
  // {
  //   testExtremeInitialValue<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolver<Group1CudaSparseMatrix>>(1, initial_value);
  //   // testExtremeInitialValue<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolver<Group20CudaSparseMatrix>>(20, initial_value);
  //   // testExtremeInitialValue<Group300CudaDenseMatrix, Group300CudaSparseMatrix, micm::CudaLinearSolver<Group300CudaSparseMatrix>>(
  //   //     300, initial_value);
  //   // testExtremeInitialValue<Group4000CudaDenseMatrix, Group4000CudaSparseMatrix, micm::CudaLinearSolver<Group4000CudaSparseMatrix>>(
  //   //     4000, initial_value);
  // }
}