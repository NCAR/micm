#include <gtest/gtest.h>

#include <micm/util/sparse_matrix.hpp>

template<class OrderingPolicy>
micm::SparseMatrix<double, OrderingPolicy> testZeroMatrix()
{
  auto builder = micm::SparseMatrix<double, OrderingPolicy>::create(3);
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 0);
  EXPECT_EQ(row_ids.size(), 0);
  EXPECT_EQ(row_starts.size(), 4);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 0);
  EXPECT_EQ(row_starts[2], 0);
  EXPECT_EQ(row_starts[3], 0);

  micm::SparseMatrix<double, OrderingPolicy> matrix{ builder };

  EXPECT_EQ(matrix.FlatBlockSize(), 0);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 0); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 0); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 3); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 3); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 0); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(1, 3); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 3); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2.0; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[1][0][0] = 2.0; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][3][0] = 2.0; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][1][1] = 2.0; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  return matrix;
}
