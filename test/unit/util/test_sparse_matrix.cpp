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
  // 0 1 0 1
  // 0 0 1 0
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
  // 0 1 0 1
  // 0 0 1 0
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
}