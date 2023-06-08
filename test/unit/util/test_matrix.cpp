#include <gtest/gtest.h>

#include <micm/util/matrix.hpp>
#include "test_matrix_policy.hpp"

TEST(Matrix, SmallMatrix)
{
  micm::Matrix<double> matrix{3, 5};

  testSmallMatrix(matrix);

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, SmallConstMatrix)
{
  micm::Matrix<double> orig_matrix{3, 5};

  testSmallConstMatrix(orig_matrix);

  const micm::Matrix<double> matrix = orig_matrix;
  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, InitializeMatrix)
{
  micm::Matrix<double> matrix{2, 3, 12.4};

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(Matrix, InitializeConstMatrix)
{
  const micm::Matrix<double> matrix{2, 3, 12.4};

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(Matrix, LoopOverMatrix)
{
  micm::Matrix<int> matrix{3, 4, 0};

  testLoopOverMatrix(matrix);
}

TEST(Matrix, LoopOverConstMatrix)
{
  micm::Matrix<int> matrix{3, 4, 0};

  testLoopOverConstMatrix(matrix);
}

TEST(Matrix, IterateOverMatrix)
{
  micm::Matrix<int> matrix{3, 4, 0};
  int i = 42;

  for (auto& elem : matrix[1]) {
    elem = i++;
  }

  const micm::Matrix<int> copy = matrix;
  
  i = 42;
  for (auto& elem : matrix[1]) {
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
  micm::Matrix<double> matrix{2, 3, 0.0};

  testConversionToVector(matrix);
}

TEST(Matrix, ConstConversionToVector)
{
  micm::Matrix<double> matrix{2, 3, 0.0};

  testConstConversionToVector(matrix);
}

TEST(Matrix, ConversionFromVector)
{
  testConversionFromVector<micm::Matrix<double>>();
}

TEST(Matrix, AssignmentFromVector)
{
  micm::Matrix matrix{4, 3, 0.0};

  testAssignmentFromVector(matrix);
}