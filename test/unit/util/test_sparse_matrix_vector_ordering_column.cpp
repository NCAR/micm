#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedColumnMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[12], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[16], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedColumnMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedColumnMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 8);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[8], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 14);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[14], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}
