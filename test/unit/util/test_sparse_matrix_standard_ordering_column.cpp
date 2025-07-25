#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

using StandardOrdering = micm::SparseMatrixStandardOrderingCompressedSparseColumn;

TEST(SparseCompressedColumnMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SingleBlockMatrix)
{
  {
    auto matrix = testSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 42;
      EXPECT_EQ(matrix.AsVector()[3], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
  }
}

TEST(SparseCompressedColumnMatrix, ConstSingleBlockMatrix)
{
  {
    auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 3);
      EXPECT_EQ(matrix.AsVector()[3], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 4);
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
  }
}

TEST(SparseCompressedColumnMatrix, MultiBlockMatrix)
{
  {
    auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(0, 2, 3);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 2, 1);
      EXPECT_EQ(elem, 12);
      matrix.AsVector()[elem] = 31;
      EXPECT_EQ(matrix.AsVector()[12], 31);
    }
  }
}

TEST(SparseCompressedColumnMatrix, Print)
{
  auto matrix = testPrint<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, PrintNonZero)
{
  testPrintNonZero<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(3, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(2, 4); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(6, 7); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
}
