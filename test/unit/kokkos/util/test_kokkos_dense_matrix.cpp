#include <gtest/gtest.h>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>

TEST(KokkosDenseMatrix, DefaultConstructor)
{
  micm::KokkosDenseMatrix<double> matrix;
  EXPECT_EQ(matrix.NumRows(), 0);
  EXPECT_EQ(matrix.NumColumns(), 0);
}

TEST(KokkosDenseMatrix, DimensionsConstructor)
{
  micm::KokkosDenseMatrix<double> matrix(3, 4);
  EXPECT_EQ(matrix.NumRows(), 3);
  EXPECT_EQ(matrix.NumColumns(), 4);
  // Padding occurs for VectorMatrix with L=4, 3 rows becomes 1 group of 4
  EXPECT_EQ(matrix.AsVector().size(), 16);
}

TEST(KokkosDenseMatrix, CopyToDeviceAndHost)
{
  micm::KokkosDenseMatrix<double> matrix(2, 2);
  matrix[0][0] = 1.0;
  matrix[0][1] = 2.0;
  matrix[1][0] = 3.0;
  matrix[1][1] = 4.0;

  matrix.CopyToDevice();

  // Clear host data manually (not using matrix.Fill(0.0) as it clears device data)
  for (auto& elem : matrix.AsVector()) elem = 0.0;
  EXPECT_EQ(matrix[0][0], 0.0);

  matrix.CopyToHost();
  EXPECT_EQ(matrix[0][0], 1.0);
  EXPECT_EQ(matrix[0][1], 2.0);
  EXPECT_EQ(matrix[1][0], 3.0);
  EXPECT_EQ(matrix[1][1], 4.0);
}
