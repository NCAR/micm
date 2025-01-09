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

TEST(LuDecompositionDoolittleInPlace, DenseMatrixStandardOrdering)
{
  testDenseMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(1);
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueStandardOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5, value);
  }
}

TEST(LuDecompositionDoolittleInPlace, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueVectorOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
  }
}