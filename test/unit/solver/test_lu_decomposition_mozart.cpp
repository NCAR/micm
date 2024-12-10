#include "test_lu_decomposition_policy.hpp"

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<double>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(LuDecompositionMozart, DenseMatrixStandardOrdering)
{
  testDenseMatrix<SparseMatrixTest, micm::LuDecompositionMozart>();
}

TEST(LuDecompositionMozart, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(1);
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionMozart>(5, value);
  }
}

TEST(LuDecompositionMozart, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>();
}

TEST(LuDecompositionMozart, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, VectorOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
  }
}