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
using Group3CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 3>;
using Group27CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 27>;
using Group32CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 32>;
using Group43CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 43>;
using Group77CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 77>;
using Group113CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 113>;
using Group193CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 193>;
using Group281CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 281>;
using Group472CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 472>;
using Group512CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 512>;
using Group739CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 739>;
using Group1130CudaDenseMatrix = micm::CudaDenseMatrix<FloatingPointType, 1130>;

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group3CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<3>>;
using Group27CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<27>>;
using Group32CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<32>>;
using Group43CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<43>>;
using Group77CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<77>>;
using Group113CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<113>>;
using Group193CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<193>>;
using Group281CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<281>>;
using Group472CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<472>>;
using Group512CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<512>>;
using Group739CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<739>>;
using Group1130CudaSparseMatrix = micm::CudaSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1130>>;

TEST(CudaLinearSolverInPlace, DenseMatrixVectorOrderingPolicy)
{
  testDenseMatrix<
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group1130CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>>();
}

TEST(CudaLinearSolverInPlace, RandomMatrixVectorOrderingPolicy)
{
  testRandomMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(
      400);
  testRandomMatrix<Group3CudaDenseMatrix, Group3CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group3CudaSparseMatrix>>(
      400);
  testRandomMatrix<Group27CudaDenseMatrix, Group27CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group27CudaSparseMatrix>>(
      400);
  testRandomMatrix<Group32CudaDenseMatrix, Group32CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group32CudaSparseMatrix>>(
      400);
  testRandomMatrix<Group43CudaDenseMatrix, Group43CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group43CudaSparseMatrix>>(
      400);
  testRandomMatrix<Group77CudaDenseMatrix, Group77CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group77CudaSparseMatrix>>(
      400);
  testRandomMatrix<
      Group113CudaDenseMatrix,
      Group113CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group113CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group193CudaDenseMatrix,
      Group193CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group193CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group281CudaDenseMatrix,
      Group281CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group281CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group472CudaDenseMatrix,
      Group472CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group472CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group512CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group739CudaDenseMatrix,
      Group739CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group739CudaSparseMatrix>>(400);
  testRandomMatrix<
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group1130CudaSparseMatrix>>(400);
}

TEST(CudaLinearSolverInPlace, DiagonalMatrixVectorOrderingPolicy)
{
  testDiagonalMatrix<Group1CudaDenseMatrix, Group1CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(
      400);
  testDiagonalMatrix<Group3CudaDenseMatrix, Group3CudaSparseMatrix, micm::CudaLinearSolverInPlace<Group3CudaSparseMatrix>>(
      400);
  testDiagonalMatrix<
      Group27CudaDenseMatrix,
      Group27CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group27CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group32CudaDenseMatrix,
      Group32CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group32CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group43CudaDenseMatrix,
      Group43CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group43CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group77CudaDenseMatrix,
      Group77CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group77CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group113CudaDenseMatrix,
      Group113CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group113CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group193CudaDenseMatrix,
      Group193CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group193CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group281CudaDenseMatrix,
      Group281CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group281CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group472CudaDenseMatrix,
      Group472CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group472CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group512CudaDenseMatrix,
      Group512CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group512CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group739CudaDenseMatrix,
      Group739CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group739CudaSparseMatrix>>(400);
  testDiagonalMatrix<
      Group1130CudaDenseMatrix,
      Group1130CudaSparseMatrix,
      micm::CudaLinearSolverInPlace<Group1130CudaSparseMatrix>>(400);
}

TEST(CudaLinearSolverInPlace, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto initial_value : initial_values)
  {
    testExtremeInitialValue<
        Group1CudaDenseMatrix,
        Group1CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group1CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group3CudaDenseMatrix,
        Group3CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group3CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group27CudaDenseMatrix,
        Group27CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group27CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group32CudaDenseMatrix,
        Group32CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group32CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group43CudaDenseMatrix,
        Group43CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group43CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group77CudaDenseMatrix,
        Group77CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group77CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group113CudaDenseMatrix,
        Group113CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group113CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group193CudaDenseMatrix,
        Group193CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group193CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group281CudaDenseMatrix,
        Group281CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group281CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group472CudaDenseMatrix,
        Group472CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group472CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group512CudaDenseMatrix,
        Group512CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group512CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group739CudaDenseMatrix,
        Group739CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group739CudaSparseMatrix>>(400, initial_value);
    testExtremeInitialValue<
        Group1130CudaDenseMatrix,
        Group1130CudaSparseMatrix,
        micm::CudaLinearSolverInPlace<Group1130CudaSparseMatrix>>(400, initial_value);
  }
}