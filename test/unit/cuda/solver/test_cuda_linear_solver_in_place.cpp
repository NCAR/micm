#include "../../solver/test_linear_solver_in_place_policy.hpp"

#include <micm/cuda/solver/cuda_linear_solver_in_place.hpp>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/linear_solver_in_place.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <random>

using FloatingPointType = double;

using Group1CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1>;
using Group20CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 20>;
using Group300CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 300>;
using Group4000CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 4000>;
using Group10000CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 10000>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group20CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<20>>;
using Group300CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<300>>;
using Group4000CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<4000>>;
using Group10000CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<10000>>;

TEST(CudaLinearSolverInPlace, DenseMatrixVectorOrderingPolicy)
{
  testDenseMatrix<
      Group1CudaDenseMatrix,
      Group1CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>>();
}

TEST(CudaLinearSolverInPlace, RandomMatrixVectorOrderingPolicy)
{
  testRandomMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(1);
  testRandomMatrix<Group20CudaDenseMatrix, Group20CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group20CudaSparseMatrix>>(
      20);
  testRandomMatrix<
      Group300CudaDenseMatrix,
      Group300CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group300CudaSparseMatrix>>(300);
  testRandomMatrix<
      Group4000CudaDenseMatrix,
      Group4000CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group4000CudaSparseMatrix>>(4000);
}

TEST(CudaLinearSolverInPlace, DiagonalMatrixVectorOrderingPolicy)
{
  testDiagonalMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(
      1);
  testDiagonalMatrix<
      Group20CudaDenseMatrix,
      Group20CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group20CudaSparseMatrix>>(20);
  testDiagonalMatrix<
      Group300CudaDenseMatrix,
      Group300CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group300CudaSparseMatrix>>(300);
  testDiagonalMatrix<
      Group4000CudaDenseMatrix,
      Group4000CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group4000CudaSparseMatrix>>(4000);
}

TEST(CudaLinearSolverInPlace, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto initial_value : initial_values)
  {
    testExtremeInitialValue<
        Group1CudaDenseMatrix,
        Group1CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(1, initial_value);
    testExtremeInitialValue<
        Group20CudaDenseMatrix,
        Group20CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group20CudaSparseMatrix>>(20, initial_value);
    testExtremeInitialValue<
        Group300CudaDenseMatrix,
        Group300CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group300CudaSparseMatrix>>(300, initial_value);
    testExtremeInitialValue<
        Group4000CudaDenseMatrix,
        Group4000CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group4000CudaSparseMatrix>>(4000, initial_value);
  }
}