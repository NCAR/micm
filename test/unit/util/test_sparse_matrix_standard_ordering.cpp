#include <gtest/gtest.h>

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include "test_sparse_matrix_policy.hpp"

using StandardOrdering = micm::SparseMatrixStandardOrdering;

TEST(SparseMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrix, SingleBlockMatrix)
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
    EXPECT_EQ(elem, 2);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[2], 21);
  }
}

TEST(SparseMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 3);
    EXPECT_EQ(matrix.AsVector()[3], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 2);
    EXPECT_EQ(matrix.AsVector()[2], 21);
  }
}

TEST(SparseMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 2);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[2], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 9);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[9], 31);
  }
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::create(3).with_element(3, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::create(3).with_element(2, 4); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::create(3).with_element(6, 7); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
}