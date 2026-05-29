#include "test_matrix_policy.hpp"

#include <micm/util/matrix.hpp>

#include <gtest/gtest.h>

TEST(Matrix, SmallMatrix)
{
  auto matrix = TestSmallMatrix<micm::Matrix>();

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, SmallConstMatrix)
{
  auto matrix = TestSmallConstMatrix<micm::Matrix>();

  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 3 * 5);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[3 * 5 - 1], 102.3);
  EXPECT_EQ(data[1 * 5 + 3], 64.7);
}

TEST(Matrix, InitializeMatrix)
{
  TestInializeMatrix<micm::Matrix>();
}

TEST(Matrix, InitializeConstMatrix)
{
  TestInializeConstMatrix<micm::Matrix>();
}

TEST(Matrix, LoopOverMatrix)
{
  TestLoopOverMatrix<micm::Matrix>();
}

TEST(Matrix, LoopOverConstMatrix)
{
  TestLoopOverConstMatrix<micm::Matrix>();
}

TEST(Matrix, IterateOverMatrix)
{
  micm::Matrix<int> matrix{ 3, 4, 0 };
  int i = 42;

  for (auto& elem : matrix[1])
  {
    elem = i++;
  }

  const micm::Matrix<int> COPY = matrix;

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

TEST(Matrix, Strides)
{
  auto matrix = TestStrides<micm::Matrix>();
  EXPECT_EQ(matrix.RowStride(), 4);
  EXPECT_EQ(matrix.ColumnStride(), 1);
}

TEST(Matrix, ConversionToVector)
{
  TestConversionToVector<micm::Matrix>();
}

TEST(Matrix, ConstConversionToVector)
{
  TestConstConversionToVector<micm::Matrix>();
}

TEST(Matrix, ConversionFromVector)
{
  TestConversionFromVector<micm::Matrix>();
}

TEST(Matrix, AssignmentFromVector)
{
  TestAssignmentFromVector<micm::Matrix>();
}

TEST(Matrix, Axpy)
{
  TestAxpy<micm::Matrix>();
}

TEST(Matrix, ForEach)
{
  TestForEach<micm::Matrix>();
}

TEST(Matrix, SetScaler)
{
  TestSetScalar<micm::Matrix>();
}

TEST(Matrix, Max)
{
  TestMax<micm::Matrix>();
}

TEST(Matrix, Min)
{
  TestMin<micm::Matrix>();
}
TEST(Matrix, ArrayFunction)
{
  TestArrayFunction<micm::Matrix>();
}

TEST(Matrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::Matrix>();
}

TEST(Matrix, MismatchedRowDimensions)
{
  TestMismatchedRowDimensions<micm::Matrix>();
}

TEST(Matrix, MismatchedColumnDimensions)
{
  TestMismatchedColumnDimensions<micm::Matrix>();
}

TEST(Matrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<micm::Matrix>();
}

TEST(Matrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::Matrix>();
}

TEST(Matrix, ColumnViewReuse)
{
  TestColumnViewReuse<micm::Matrix>();
}

TEST(Matrix, FunctionReusability)
{
  TestFunctionReusability<micm::Matrix>();
}

TEST(Matrix, ConstMatrixFunction)
{
  TestConstMatrixFunction<micm::Matrix>();
}

TEST(Matrix, EmptyMatrixFunction)
{
  TestEmptyMatrixFunction<micm::Matrix>();
}

// Flexible row count Tests
TEST(Matrix, MultiMatrixDifferentRowsFromCreation)
{
  TestMultiMatrixDifferentRowsFromCreation<micm::Matrix>();
}

TEST(Matrix, MatrixVectorDifferentRowsFromCreation)
{
  TestMatrixVectorDifferentRowsFromCreation<micm::Matrix>();
}

TEST(Matrix, MismatchedRowsAtInvocation)
{
  TestMismatchedRowsAtInvocation<micm::Matrix>();
}

TEST(Matrix, MultipleMatricesMismatchedRowsAtInvocation)
{
  TestMultipleMatricesMismatchedRowsAtInvocation<micm::Matrix>();
}

TEST(Matrix, WrongColumnCountAtInvocation)
{
  TestWrongColumnCountAtInvocation<micm::Matrix>();
}

TEST(Matrix, Print)
{
  TestPrint<micm::Matrix>();
}
// Vector support Tests
TEST(Matrix, VectorInMatrixFunction)
{
  TestVectorInMatrixFunction<micm::Matrix>();
}

TEST(Matrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::Matrix>();
}

TEST(Matrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::Matrix>();
}

TEST(Matrix, EmptyVectorNonEmptyMatrix)
{
  TestEmptyVectorNonEmptyMatrix<micm::Matrix>();
}

TEST(Matrix, NonEmptyVectorEmptyMatrix)
{
  TestNonEmptyVectorEmptyMatrix<micm::Matrix>();
}

TEST(Matrix, EmptyVectorEmptyMatrix)
{
  TestEmptyVectorEmptyMatrix<micm::Matrix>();
}

TEST(Matrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::Matrix>();
}

TEST(Matrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::Matrix>();
}

TEST(Matrix, MultipleMatricesOneVector)
{
  TestMultipleMatricesOneVector<micm::Matrix>();
}

TEST(Matrix, MultipleMatricesDifferentRowsVector)
{
  TestMultipleMatricesDifferentRowsVector<micm::Matrix>();
}

TEST(Matrix, VectorSizeMatchesOneMatrixOnly)
{
  TestVectorSizeMatchesOneMatrixOnly<micm::Matrix>();
}

TEST(Matrix, ConstVector)
{
  TestConstVector<micm::Matrix>();
}

TEST(Matrix, MutableVector)
{
  TestMutableVector<micm::Matrix>();
}

TEST(Matrix, FunctionReusabilityWithVectors)
{
  TestFunctionReusabilityWithVectors<micm::Matrix>();
}

TEST(Matrix, FunctionInvocationWithWrongSizedVector)
{
  TestFunctionInvocationWithWrongSizedVector<micm::Matrix>();
}

TEST(Matrix, ArraySupport)
{
  TestArraySupport<micm::Matrix>();
}

TEST(Matrix, MixedVectorColumnViewRowVariable)
{
  TestMixedVectorColumnViewRowVariable<micm::Matrix>();
}

TEST(Matrix, IntegerVector)
{
  TestIntegerVector<micm::Matrix>();
}

TEST(Matrix, FunctionWithConstSignature)
{
  TestFunctionWithConstSignature<micm::Matrix>();
}
