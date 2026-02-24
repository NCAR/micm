#include "test_matrix_policy.hpp"

#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<class T>
using Group1MatrixAlias = micm::VectorMatrix<T, 1>;
template<class T>
using Group2MatrixAlias = micm::VectorMatrix<T, 2>;
template<class T>
using Group3MatrixAlias = micm::VectorMatrix<T, 3>;
template<class T>
using Group4MatrixAlias = micm::VectorMatrix<T, 4>;

TEST(VectorMatrix, SmallVectorMatrix)
{
  auto matrix = testSmallMatrix<Group2MatrixAlias>();

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 2);
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 * 5 + 0 + 2 * 4], 102.3);
  EXPECT_EQ(data[1 + 2 * 3], 64.7);
}

TEST(VectorMatrix, SmallConstVectorMatrix)
{
  auto matrix = testSmallConstMatrix<Group4MatrixAlias>();

  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 4 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 1);
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 + 4 * 4], 102.3);
  EXPECT_EQ(data[1 + 4 * 3], 64.7);
}

TEST(VectorMatrix, InitializeVectorMatrix)
{
  testInializeMatrix<Group1MatrixAlias>();
}

TEST(VectorMatrix, InitializeConstVectorMatrix)
{
  testInializeConstMatrix<Group2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverVectorMatrix)
{
  testLoopOverMatrix<Group2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverConstVectorMatrix)
{
  testLoopOverConstMatrix<Group1MatrixAlias>();
}

TEST(VectorMatrix, Strides)
{
  auto matrix3vec = testStrides<Group3MatrixAlias>();
  EXPECT_EQ(matrix3vec.RowStride(), 1);
  EXPECT_EQ(matrix3vec.ColumnStride(), 3);
  auto matrix4vec = testStrides<Group4MatrixAlias>();
  EXPECT_EQ(matrix4vec.RowStride(), 1);
  EXPECT_EQ(matrix4vec.ColumnStride(), 4);
}

TEST(VectorMatrix, ConversionToVector)
{
  testConversionToVector<Group3MatrixAlias>();
}

TEST(VectorMatrix, ConstConversionToVector)
{
  testConstConversionToVector<Group1MatrixAlias>();
}

TEST(VectorMatrix, ConversionFromVector)
{
  testConversionFromVector<Group2MatrixAlias>();
}

TEST(VectorMatrix, AssignmentFromVector)
{
  testAssignmentFromVector<Group2MatrixAlias>();
}

TEST(VectorMatrix, Axpy)
{
  testAxpy<Group1MatrixAlias>();
  testAxpy<Group2MatrixAlias>();
  testAxpy<Group3MatrixAlias>();
  testAxpy<Group4MatrixAlias>();
}

TEST(VectorMatrix, ForEach)
{
  testForEach<Group1MatrixAlias>();
  testForEach<Group2MatrixAlias>();
  testForEach<Group3MatrixAlias>();
  testForEach<Group4MatrixAlias>();
}

TEST(VectorMatrix, SetScaler)
{
  testSetScalar<Group1MatrixAlias>();
  testSetScalar<Group2MatrixAlias>();
  testSetScalar<Group3MatrixAlias>();
  testSetScalar<Group4MatrixAlias>();
}

TEST(VectorMatrix, Max)
{
  testMax<Group1MatrixAlias>();
  testMax<Group2MatrixAlias>();
  testMax<Group3MatrixAlias>();
  testMax<Group4MatrixAlias>();
}

TEST(VectorMatrix, Min)
{
  testMin<Group1MatrixAlias>();
  testMin<Group2MatrixAlias>();
  testMin<Group3MatrixAlias>();
  testMin<Group4MatrixAlias>();
}

TEST(VectorMatrix, Print)
{
  testPrint<Group1MatrixAlias>();
  testPrint<Group2MatrixAlias>();
  testPrint<Group3MatrixAlias>();
  testPrint<Group4MatrixAlias>();
}

TEST(VectorMatrix, ArrayFunction)
{
  testArrayFunction<Group1MatrixAlias>();
  testArrayFunction<Group2MatrixAlias>();
  testArrayFunction<Group3MatrixAlias>();
  testArrayFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<Group1MatrixAlias>();
  testMultiMatrixArrayFunction<Group2MatrixAlias>();
  testMultiMatrixArrayFunction<Group3MatrixAlias>();
  testMultiMatrixArrayFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, MismatchedRowDimensions)
{
  testMismatchedRowDimensions<Group1MatrixAlias>();
  testMismatchedRowDimensions<Group2MatrixAlias>();
  testMismatchedRowDimensions<Group3MatrixAlias>();
  testMismatchedRowDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, MismatchedColumnDimensions)
{
  testMismatchedColumnDimensions<Group1MatrixAlias>();
  testMismatchedColumnDimensions<Group2MatrixAlias>();
  testMismatchedColumnDimensions<Group3MatrixAlias>();
  testMismatchedColumnDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<Group1MatrixAlias>();
  testWrongMatrixDimensions<Group2MatrixAlias>();
  testWrongMatrixDimensions<Group3MatrixAlias>();
  testWrongMatrixDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<Group1MatrixAlias>();
  testMultipleTemporaries<Group2MatrixAlias>();
  testMultipleTemporaries<Group3MatrixAlias>();
  testMultipleTemporaries<Group4MatrixAlias>();
}

TEST(VectorMatrix, ColumnViewReuse)
{
  testColumnViewReuse<Group1MatrixAlias>();
  testColumnViewReuse<Group2MatrixAlias>();
  testColumnViewReuse<Group3MatrixAlias>();
  testColumnViewReuse<Group4MatrixAlias>();
}

TEST(VectorMatrix, FunctionReusability)
{
  testFunctionReusability<Group1MatrixAlias>();
  testFunctionReusability<Group2MatrixAlias>();
  testFunctionReusability<Group3MatrixAlias>();
  testFunctionReusability<Group4MatrixAlias>();
}

TEST(VectorMatrix, ConstMatrixFunction)
{
  testConstMatrixFunction<Group1MatrixAlias>();
  testConstMatrixFunction<Group2MatrixAlias>();
  testConstMatrixFunction<Group3MatrixAlias>();
  testConstMatrixFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, EmptyMatrixFunction)
{
  testEmptyMatrixFunction<Group1MatrixAlias>();
  testEmptyMatrixFunction<Group2MatrixAlias>();
  testEmptyMatrixFunction<Group3MatrixAlias>();
  testEmptyMatrixFunction<Group4MatrixAlias>();
}
