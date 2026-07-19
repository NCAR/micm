#include "test_lu_decomposition_in_place_policy.hpp"

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<micm::Real>;

using Group1SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<4>>;

TEST(LuDecompositionMozartInPlace, DenseMatrixStandardOrdering)
{
  TestDenseMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>();
}

TEST(LuDecompositionMozartInPlace, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(1);
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, AgnosticToInitialValueStandardOrdering)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionMozartInPlace>(5, value);
  }
}

TEST(LuDecompositionMozartInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  TestDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  TestDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
  TestDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>();
}

TEST(LuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
  TestDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5);
}

TEST(LuDecompositionMozartInPlace, AgnosticToInitialValueVectorOrdering)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    TestExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    TestExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
    TestExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionMozartInPlace>(5, value);
  }
}