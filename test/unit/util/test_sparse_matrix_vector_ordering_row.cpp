#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedRowMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();

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
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

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
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

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
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}

TEST(SparseVectorCompressedRowMatrix, Print)
{
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}
