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

TEST(LuDecomposition, DenseMatrixStandardOrdering)
{
  testDenseMatrix<SparseMatrixTest, micm::LuDecomposition>();
}

TEST(LuDecomposition, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecomposition>(1);
  testRandomMatrix<SparseMatrixTest, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<SparseMatrixTest, micm::LuDecomposition>(5, value);
  }
}

TEST(LuDecomposition, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>();
}

TEST(LuDecomposition, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, VectorOrderingAgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::LuDecomposition>(5, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::LuDecomposition>(5, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::LuDecomposition>(5, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::LuDecomposition>(5, value);
  }
}