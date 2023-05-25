#include <gtest/gtest.h>

#include <micm/util/matrix.hpp>

TEST(Matrix, SmallMatrix)
{
  micm::Matrix<double> matrix{3, 5};

  matrix[1][3] = 64.7;
  matrix[0][0] = 41.2;
  matrix[2][4] = 102.3;

  EXPECT_EQ(matrix[1][3], 64.7);
  EXPECT_EQ(matrix[0][0], 41.2);
  EXPECT_EQ(matrix[2][4], 102.3);

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}