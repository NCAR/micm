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
    std::size_t elem = matrix_csr.VectorIndex(3, 2);
    EXPECT_EQ(elem, 16);
    matrix_csr.AsVector()[elem] = 42;
    EXPECT_EQ(matrix_csr.AsVector()[16], 42);
  }
  {
    std::size_t elem = matrix_csr.VectorIndex(2, 3);
    EXPECT_EQ(elem, 12);
    matrix_csr.AsVector()[elem] = 21;
    EXPECT_EQ(matrix_csr.AsVector()[12], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    std::size_t elem = matrix_csr.VectorIndex(3, 2);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix_csr.AsVector()[8], 42);
  }
  {
    std::size_t elem = matrix_csr.VectorIndex(2, 3);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix_csr.AsVector()[6], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    std::size_t elem = matrix_csr.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 6);
    matrix_csr.AsVector()[elem] = 21;
    EXPECT_EQ(matrix_csr.AsVector()[6], 21);
  }
  {
    std::size_t elem = matrix_csr.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 14);
    matrix_csr.AsVector()[elem] = 31;
    EXPECT_EQ(matrix_csr.AsVector()[14], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}
