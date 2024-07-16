#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<2>>();
}

TEST(SparseVectorMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(SparseVectorMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(SparseVectorMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<4>>();
}

TEST(SparseVectorMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[16], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[12], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 6);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 14);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[14], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}
