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
  testDenseMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>();
}

TEST(LuDecompositionDoolittle, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(1);
  testRandomMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<SparseMatrixTest, micm::LuDecompositionDoolittle>(5, value);
  }
}

TEST(LuDecompositionDoolittle, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>();
}

TEST(LuDecompositionDoolittle, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5);
}

TEST(LuDecompositionDoolittle, VectorOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecompositionDoolittle>(5, value);
  }
}