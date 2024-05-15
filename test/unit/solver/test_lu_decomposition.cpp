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

TEST(LuDecomposition, SingularMatrixStandardOrdering)
{
  testSingularMatrix<SparseMatrixTest, micm::LuDecomposition>();
}

TEST(LuDecomposition, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecomposition>(5);
}

TEST(LuDecomposition, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>();
}

TEST(LuDecomposition, SingluarMatrixVectorOrdering)
{
  testSingularMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>();
  testSingularMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>();
  testSingularMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>();
  testSingularMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>();
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
