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
  TestDenseMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(1);
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueStandardOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionDoolittleInPlace>(5, value);
  }
}

TEST(LuDecompositionDoolittleInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  TestDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  TestDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
  TestDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
  TestDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueVectorOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    TestExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    TestExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
    TestExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionDoolittleInPlace>(5, value);
  }
}