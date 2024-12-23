#include "test_lu_decomposition_in_place_policy.hpp"

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<double>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(LuDecompositionMozartInPlace, DenseMatrixStandardOrdering)
{
  testDenseMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>();
}

TEST(LuDecompositionMozartInPlace, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(1);
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, AgnosticToInitialValueStandardOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5, value);
  }
}

TEST(LuDecompositionMozartInPlace, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
}

TEST(LuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, AgnosticToInitialValueVectorOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
  }
}