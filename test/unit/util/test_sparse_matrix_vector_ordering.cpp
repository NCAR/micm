#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorMatrix, SingleBlockMatrix)
{
  auto matrix_csr = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
  auto matrix_csc = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();

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
  EXPECT_EQ(matrix_csr.GroupVectorSize(), 4);
  EXPECT_EQ(matrix_csr.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix_csr.NumberOfGroups(1), 1);
  EXPECT_EQ(matrix_csc.GroupVectorSize(), 4);
  EXPECT_EQ(matrix_csc.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix_csc.NumberOfGroups(1), 1);
}

TEST(SparseVectorMatrix, ConstSingleBlockMatrix)
{
  auto matrix_csr = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  auto matrix_csc = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

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
  EXPECT_EQ(matrix_csr.GroupVectorSize(), 2);
  EXPECT_EQ(matrix_csr.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix_csr.NumberOfGroups(1), 1);
  EXPECT_EQ(matrix_csc.GroupVectorSize(), 2);
  EXPECT_EQ(matrix_csc.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix_csc.NumberOfGroups(1), 1);
}

TEST(SparseVectorMatrix, MultiBlockMatrix)
{
  auto matrix_csr = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  auto matrix_csc = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

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
  EXPECT_EQ(matrix_csr.GroupVectorSize(), 2);
  EXPECT_EQ(matrix_csr.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix_csr.NumberOfGroups(4), 2);
  EXPECT_EQ(matrix_csc.GroupVectorSize(), 2);
  EXPECT_EQ(matrix_csc.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix_csc.NumberOfGroups(4), 2);
}
