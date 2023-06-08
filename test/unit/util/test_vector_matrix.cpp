#include <gtest/gtest.h>

#include <micm/util/vector_matrix.hpp>
#include "test_matrix_policy.hpp"

template<class T>
using Block1MatrixAlias = micm::VectorMatrix<T, 1>;
template<class T>
using Block2MatrixAlias = micm::VectorMatrix<T, 2>;
template<class T>
using Block3MatrixAlias = micm::VectorMatrix<T, 3>;
template<class T>
using Block4MatrixAlias = micm::VectorMatrix<T, 4>;

TEST(VectorMatrix, SmallVectorMatrix)
{
  auto matrix = testSmallMatrix<Block2MatrixAlias>();

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 * 5 + 0 + 2 * 4], 102.3);
  EXPECT_EQ(data[1 + 2 * 3], 64.7);
}

TEST(VectorMatrix, SmallConstVectorMatrix)
{
  auto matrix = testSmallConstMatrix<Block4MatrixAlias>();

  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 + 4 * 4], 102.3);
  EXPECT_EQ(data[1 + 4 * 3], 64.7);
}

TEST(VectorMatrix, InitializeVectorMatrix)
{
  testInializeMatrix<Block1MatrixAlias>();
}

TEST(VectorMatrix, InitializeConstVectorMatrix)
{
  testInializeConstMatrix<Block2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverVectorMatrix)
{
  testLoopOverMatrix<Block2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverConstVectorMatrix)
{
  testLoopOverConstMatrix<Block1MatrixAlias>();
}

TEST(VectorMatrix, ConversionToVector)
{
  testConversionToVector<Block3MatrixAlias>();
}

TEST(VectorMatrix, ConstConversionToVector)
{
  testConstConversionToVector<Block1MatrixAlias>();
}

TEST(VectorMatrix, ConversionFromVector)
{
  testConversionFromVector<Block2MatrixAlias>();
}

TEST(VectorMatrix, AssignmentFromVector)
{
  testAssignmentFromVector<Block2MatrixAlias>();
}