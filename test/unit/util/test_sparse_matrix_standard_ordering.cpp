#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

using StandardOrderingCSR = micm::SparseMatrixStandardOrderingCompressedSparseRow;
using StandardOrderingCSC = micm::SparseMatrixStandardOrderingCompressedSparseColumn;

TEST(SparseMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, StandardOrderingCSR>();
  testZeroMatrix<micm::SparseMatrix, StandardOrderingCSC>();
}

TEST(SparseMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, StandardOrderingCSR>();
  testConstZeroMatrix<micm::SparseMatrix, StandardOrderingCSC>();
}

TEST(SparseMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, StandardOrderingCSR>();
  testSetScalar<micm::SparseMatrix, StandardOrderingCSC>();
}

TEST(SparseMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, StandardOrderingCSR>();
  testAddToDiagonal<micm::SparseMatrix, StandardOrderingCSC>();
}

TEST(SparseMatrix, SingleBlockMatrix)
{
  {
    auto matrix = testSingleBlockMatrix<micm::SparseMatrix, StandardOrderingCSR>();
  
    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 42;
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
  {
    auto matrix = testSingleBlockMatrix<micm::SparseMatrix, StandardOrderingCSC>();
  
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

TEST(SparseMatrix, ConstSingleBlockMatrix)
{
  {
    auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrderingCSR>();
    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
  {
    auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrderingCSC>();
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

TEST(SparseMatrix, MultiBlockMatrix)
{
  {
    auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrderingCSR>();

    {
      std::size_t elem = matrix.VectorIndex(0, 2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 2, 1);
      EXPECT_EQ(elem, 12);
      matrix.AsVector()[elem] = 31;
      EXPECT_EQ(matrix.AsVector()[12], 31);
    }
  }
  {
    auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrderingCSC>();

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