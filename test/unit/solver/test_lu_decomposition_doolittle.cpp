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

TEST(LuDecompositionDoolittle, DenseMatrixStandardOrdering)
{
  TestDenseMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>();
}

TEST(LuDecompositionDoolittle, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(1);
  TestRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, AgnosticToInitialValueStandardOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionDoolittle>(5, value);
  }
}

TEST(LuDecompositionDoolittle, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  TestDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  TestDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  TestDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>();
}

TEST(LuDecompositionDoolittle, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  TestDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, AgnosticToInitialValueVectorOrdering)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    TestExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    TestExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    TestExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
  }
}