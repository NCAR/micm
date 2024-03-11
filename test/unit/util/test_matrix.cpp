#include <gtest/gtest.h>

#include <micm/util/matrix.hpp>

#include "test_matrix_policy.hpp"

TEST(Matrix, SmallMatrix)
{
  auto matrix = testSmallMatrix<micm::Matrix>();

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, SmallConstMatrix)
{
  auto matrix = testSmallConstMatrix<micm::Matrix>();

  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, InitializeMatrix)
{
  testInializeMatrix<micm::Matrix>();
}

TEST(Matrix, InitializeConstMatrix)
{
  testInializeConstMatrix<micm::Matrix>();
}

TEST(Matrix, LoopOverMatrix)
{
  testLoopOverMatrix<micm::Matrix>();
}

TEST(Matrix, LoopOverConstMatrix)
{
  testLoopOverConstMatrix<micm::Matrix>();
}

TEST(Matrix, IterateOverMatrix)
{
  micm::Matrix<int> matrix{ 3, 4, 0 };
  int i = 42;

  for (auto& elem : matrix[1])
  {
    elem = i++;
  }

  const micm::Matrix<int> copy = matrix;

  i = 42;
  for (auto& elem : matrix[1])
  {
    EXPECT_EQ(elem, i++);
  }

  EXPECT_EQ(matrix[0][3], 0);
  EXPECT_EQ(matrix[1][0], 42);
  EXPECT_EQ(matrix[1][1], 43);
  EXPECT_EQ(matrix[1][3], 45);
  EXPECT_EQ(matrix[2][0], 0);
}

TEST(Matrix, ConversionToVector)
{
  testConversionToVector<micm::Matrix>();
}

TEST(Matrix, ConstConversionToVector)
{
  testConstConversionToVector<micm::Matrix>();
}

TEST(Matrix, ConversionFromVector)
{
  testConversionFromVector<micm::Matrix>();
}

TEST(Matrix, AssignmentFromVector)
{
  testAssignmentFromVector<micm::Matrix>();
}

TEST(Matrix, ForEach)
{
  testForEach<micm::Matrix>();
}

TEST(Matrix, SetScaler)
{
  testSetScalar<micm::Matrix>();
}