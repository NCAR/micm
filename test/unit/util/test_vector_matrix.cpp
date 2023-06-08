#include <gtest/gtest.h>

#include <micm/util/vector_matrix.hpp>
#include "test_matrix_policy.hpp"

TEST(VectorMatrix, SmallVectorMatrix)
{
  micm::VectorMatrix<double, 2> matrix{3, 5};

  testSmallMatrix(matrix);

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 * 5 + 0 + 2 * 4], 102.3);
  EXPECT_EQ(data[1 + 2 * 3], 64.7);
}

TEST(VectorMatrix, SmallConstVectorMatrix)
{
  micm::VectorMatrix<double, 4> orig_matrix{3, 5};

  testSmallConstMatrix(orig_matrix);

  const micm::VectorMatrix<double, 4> matrix = orig_matrix;
  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 + 4 * 4], 102.3);
  EXPECT_EQ(data[1 + 4 * 3], 64.7);
}

TEST(VectorMatrix, InitializeVectorMatrix)
{
  micm::VectorMatrix<double, 1> matrix{2, 3, 12.4};

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(VectorMatrix, InitializeConstVectorMatrix)
{
  const micm::VectorMatrix<double, 2> matrix{2, 3, 12.4};

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(VectorMatrix, LoopOverVectorMatrix)
{
  micm::VectorMatrix<int, 2> matrix{3, 4, 0};

  testLoopOverMatrix(matrix);
}

TEST(VectorMatrix, LoopOverConstVectorMatrix)
{
  micm::VectorMatrix<int, 1> matrix{3, 4, 0};

  testLoopOverConstMatrix(matrix);
}

TEST(VectorMatrix, ConversionToVector)
{
  micm::VectorMatrix<double, 3> matrix{2, 3, 0.0};

  testConversionToVector(matrix);
}

TEST(VectorMatrix, ConstConversionToVector)
{
  micm::VectorMatrix<double, 1> matrix{2, 3, 0.0};

  testConstConversionToVector(matrix);
}

TEST(VectorMatrix, ConversionFromVector)
{
  testConversionFromVector<micm::VectorMatrix<double, 2>>();
}

TEST(VectorMatrix, AssignmentFromVector)
{
  micm::VectorMatrix<double, 2> matrix{4, 3, 0.0};

  testAssignmentFromVector(matrix);
}