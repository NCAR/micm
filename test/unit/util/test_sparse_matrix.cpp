#include <gtest/gtest.h>

#include <micm/util/sparse_matrix.hpp>

TEST(SparseMatrix, ZeroMatrix)
{
  auto builder = micm::SparseMatrix<double>::create(3);
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 0);
  EXPECT_EQ(row_ids.size(), 0);
  EXPECT_EQ(row_starts.size(), 4);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 0);
  EXPECT_EQ(row_starts[2], 0);
  EXPECT_EQ(row_starts[3], 0);

  micm::SparseMatrix<double> matrix{ builder };

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
}

TEST(SparseMatrix, SingleBlockMatrix)
{
  auto builder = micm::SparseMatrix<int>::create(4)
                     .with_element(0, 1)
                     .with_element(3, 2)
                     .with_element(0, 1)
                     .with_element(2, 3)
                     .with_element(2, 1);
  // 0 X 0 0
  // 0 0 0 0
  // 0 X 0 X
  // 0 0 X 0
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 4);
  EXPECT_EQ(row_ids.size(), 4);
  EXPECT_EQ(row_ids[0], 1);
  EXPECT_EQ(row_ids[1], 1);
  EXPECT_EQ(row_ids[2], 3);
  EXPECT_EQ(row_ids[3], 2);
  EXPECT_EQ(row_starts.size(), 5);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 1);
  EXPECT_EQ(row_starts[2], 1);
  EXPECT_EQ(row_starts[3], 3);
  EXPECT_EQ(row_starts[4], 4);

  micm::SparseMatrix<int> matrix{ builder };

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

  matrix[0][2][1] = 45;
  EXPECT_EQ(matrix[0][2][1], 45);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(4, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 5); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 0, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 1); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[1][0][0] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][5][0] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][1][1] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
}

TEST(SparseMatrix, MultiBlockMatrix)
{
  auto builder = micm::SparseMatrix<int>::create(4)
                     .with_element(0, 1)
                     .with_element(3, 2)
                     .with_element(0, 1)
                     .with_element(2, 3)
                     .with_element(2, 1)
                     .number_of_blocks(3);
  // 0 X 0 0
  // 0 0 0 0
  // 0 X 0 X
  // 0 0 X 0
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 4 * 3);
  EXPECT_EQ(row_ids.size(), 4);
  EXPECT_EQ(row_ids[0], 1);
  EXPECT_EQ(row_ids[1], 1);
  EXPECT_EQ(row_ids[2], 3);
  EXPECT_EQ(row_ids[3], 2);
  EXPECT_EQ(row_starts.size(), 5);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 1);
  EXPECT_EQ(row_starts[2], 1);
  EXPECT_EQ(row_starts[3], 3);
  EXPECT_EQ(row_starts[4], 4);

  micm::SparseMatrix<int> matrix{ builder };

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

  matrix[0][2][1] = 45;
  EXPECT_EQ(matrix[0][2][1], 45);

  matrix[2][2][3] = 64;
  EXPECT_EQ(matrix[2][2][3], 64);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 4, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 1, 5); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(4, 0, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 1, 1); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 0, 2); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 1); } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "Multi-block SparseMatrix access must specify block index");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[3][0][0] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][5][0] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try { matrix[0][1][1] = 2; } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix zero element access not allowed");
        throw;
      },
      std::invalid_argument);
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try {
        auto builder = micm::SparseMatrix<double>::create(3).with_element(3, 0);
      } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try {
        auto builder = micm::SparseMatrix<double>::create(3).with_element(2, 4);
      } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
  EXPECT_THROW(
      try {
        auto builder = micm::SparseMatrix<double>::create(3).with_element(6, 7);
      } catch (const std::invalid_argument& e) {
        EXPECT_STREQ(e.what(), "SparseMatrix element out of range");
        throw;
      },
      std::invalid_argument);
}