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
  TestDenseMatrix<SparseMatrixTest, micm::LuDecompositionMozart>();
}

TEST(LuDecompositionMozart, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(1);
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, AgnosticToInitialValueStandardOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionMozart>(5, value);
  }
}

TEST(LuDecompositionMozart, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>();
  TestDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>();
  TestDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>();
  TestDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>();
}

TEST(LuDecompositionMozart, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5);
  TestDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5);
}

TEST(LuDecompositionMozart, AgnosticToInitialValueVectorOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    TestExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    TestExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
    TestExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionMozart>(5, value);
  }
}