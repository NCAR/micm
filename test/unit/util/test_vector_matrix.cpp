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
  auto matrix = TestSmallMatrix<Group2MatrixAlias>();

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
  auto matrix = TestSmallConstMatrix<Group4MatrixAlias>();

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
  TestInializeMatrix<Group1MatrixAlias>();
}

TEST(VectorMatrix, InitializeConstVectorMatrix)
{
  TestInializeConstMatrix<Group2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverVectorMatrix)
{
  TestLoopOverMatrix<Group2MatrixAlias>();
}

TEST(VectorMatrix, LoopOverConstVectorMatrix)
{
  TestLoopOverConstMatrix<Group1MatrixAlias>();
}

TEST(VectorMatrix, Strides)
{
  auto matrix3vec = TestStrides<Group3MatrixAlias>();
  EXPECT_EQ(matrix3vec.RowStride(), 1);
  EXPECT_EQ(matrix3vec.ColumnStride(), 3);
  auto matrix4vec = TestStrides<Group4MatrixAlias>();
  EXPECT_EQ(matrix4vec.RowStride(), 1);
  EXPECT_EQ(matrix4vec.ColumnStride(), 4);
}

TEST(VectorMatrix, ConversionToVector)
{
  TestConversionToVector<Group3MatrixAlias>();
}

TEST(VectorMatrix, ConstConversionToVector)
{
  TestConstConversionToVector<Group1MatrixAlias>();
}

TEST(VectorMatrix, ConversionFromVector)
{
  TestConversionFromVector<Group2MatrixAlias>();
}

TEST(VectorMatrix, AssignmentFromVector)
{
  TestAssignmentFromVector<Group2MatrixAlias>();
}

TEST(VectorMatrix, Axpy)
{
  TestAxpy<Group1MatrixAlias>();
  TestAxpy<Group2MatrixAlias>();
  TestAxpy<Group3MatrixAlias>();
  TestAxpy<Group4MatrixAlias>();
}

TEST(VectorMatrix, ForEach)
{
  TestForEach<Group1MatrixAlias>();
  TestForEach<Group2MatrixAlias>();
  TestForEach<Group3MatrixAlias>();
  TestForEach<Group4MatrixAlias>();
}

TEST(VectorMatrix, SetScaler)
{
  TestSetScalar<Group1MatrixAlias>();
  TestSetScalar<Group2MatrixAlias>();
  TestSetScalar<Group3MatrixAlias>();
  TestSetScalar<Group4MatrixAlias>();
}

TEST(VectorMatrix, Max)
{
  TestMax<Group1MatrixAlias>();
  TestMax<Group2MatrixAlias>();
  TestMax<Group3MatrixAlias>();
  TestMax<Group4MatrixAlias>();
}

TEST(VectorMatrix, Min)
{
  TestMin<Group1MatrixAlias>();
  TestMin<Group2MatrixAlias>();
  TestMin<Group3MatrixAlias>();
  TestMin<Group4MatrixAlias>();
}

TEST(VectorMatrix, Print)
{
  TestPrint<Group1MatrixAlias>();
  TestPrint<Group2MatrixAlias>();
  TestPrint<Group3MatrixAlias>();
  TestPrint<Group4MatrixAlias>();
}

TEST(VectorMatrix, ArrayFunction)
{
  TestArrayFunction<Group1MatrixAlias>();
  TestArrayFunction<Group2MatrixAlias>();
  TestArrayFunction<Group3MatrixAlias>();
  TestArrayFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<Group1MatrixAlias>();
  TestMultiMatrixArrayFunction<Group2MatrixAlias>();
  TestMultiMatrixArrayFunction<Group3MatrixAlias>();
  TestMultiMatrixArrayFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, MismatchedRowDimensions)
{
  TestMismatchedRowDimensions<Group1MatrixAlias>();
  TestMismatchedRowDimensions<Group2MatrixAlias>();
  TestMismatchedRowDimensions<Group3MatrixAlias>();
  TestMismatchedRowDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, MismatchedColumnDimensions)
{
  TestMismatchedColumnDimensions<Group1MatrixAlias>();
  TestMismatchedColumnDimensions<Group2MatrixAlias>();
  TestMismatchedColumnDimensions<Group3MatrixAlias>();
  TestMismatchedColumnDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<Group1MatrixAlias>();
  TestWrongMatrixDimensions<Group2MatrixAlias>();
  TestWrongMatrixDimensions<Group3MatrixAlias>();
  TestWrongMatrixDimensions<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<Group1MatrixAlias>();
  TestMultipleTemporaries<Group2MatrixAlias>();
  TestMultipleTemporaries<Group3MatrixAlias>();
  TestMultipleTemporaries<Group4MatrixAlias>();
}

TEST(VectorMatrix, ColumnViewReuse)
{
  TestColumnViewReuse<Group1MatrixAlias>();
  TestColumnViewReuse<Group2MatrixAlias>();
  TestColumnViewReuse<Group3MatrixAlias>();
  TestColumnViewReuse<Group4MatrixAlias>();
}

TEST(VectorMatrix, FunctionReusability)
{
  TestFunctionReusability<Group1MatrixAlias>();
  TestFunctionReusability<Group2MatrixAlias>();
  TestFunctionReusability<Group3MatrixAlias>();
  TestFunctionReusability<Group4MatrixAlias>();
}

TEST(VectorMatrix, ConstMatrixFunction)
{
  TestConstMatrixFunction<Group1MatrixAlias>();
  TestConstMatrixFunction<Group2MatrixAlias>();
  TestConstMatrixFunction<Group3MatrixAlias>();
  TestConstMatrixFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, EmptyMatrixFunction)
{
  TestEmptyMatrixFunction<Group1MatrixAlias>();
  TestEmptyMatrixFunction<Group2MatrixAlias>();
  TestEmptyMatrixFunction<Group3MatrixAlias>();
  TestEmptyMatrixFunction<Group4MatrixAlias>();
}

// Flexible row count Tests
TEST(VectorMatrix, MultiMatrixDifferentRowsFromCreation)
{
  TestMultiMatrixDifferentRowsFromCreation<Group1MatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group2MatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group3MatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group4MatrixAlias>();
}

TEST(VectorMatrix, MatrixVectorDifferentRowsFromCreation)
{
  TestMatrixVectorDifferentRowsFromCreation<Group1MatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group2MatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group3MatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group4MatrixAlias>();
}

TEST(VectorMatrix, MismatchedRowsAtInvocation)
{
  TestMismatchedRowsAtInvocation<Group1MatrixAlias>();
  TestMismatchedRowsAtInvocation<Group2MatrixAlias>();
  TestMismatchedRowsAtInvocation<Group3MatrixAlias>();
  TestMismatchedRowsAtInvocation<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleMatricesMismatchedRowsAtInvocation)
{
  TestMultipleMatricesMismatchedRowsAtInvocation<Group1MatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group2MatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group3MatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group4MatrixAlias>();
}

TEST(VectorMatrix, WrongColumnCountAtInvocation)
{
  TestWrongColumnCountAtInvocation<Group1MatrixAlias>();
  TestWrongColumnCountAtInvocation<Group2MatrixAlias>();
  TestWrongColumnCountAtInvocation<Group3MatrixAlias>();
  TestWrongColumnCountAtInvocation<Group4MatrixAlias>();
}

// Vector support Tests
TEST(VectorMatrix, VectorInMatrixFunction)
{
  TestVectorInMatrixFunction<Group1MatrixAlias>();
  TestVectorInMatrixFunction<Group2MatrixAlias>();
  TestVectorInMatrixFunction<Group3MatrixAlias>();
  TestVectorInMatrixFunction<Group4MatrixAlias>();
}

TEST(VectorMatrix, VectorTooSmall)
{
  TestVectorTooSmall<Group1MatrixAlias>();
  TestVectorTooSmall<Group2MatrixAlias>();
  TestVectorTooSmall<Group3MatrixAlias>();
  TestVectorTooSmall<Group4MatrixAlias>();
}

TEST(VectorMatrix, VectorTooLarge)
{
  TestVectorTooLarge<Group1MatrixAlias>();
  TestVectorTooLarge<Group2MatrixAlias>();
  TestVectorTooLarge<Group3MatrixAlias>();
  TestVectorTooLarge<Group4MatrixAlias>();
}

TEST(VectorMatrix, EmptyVectorNonEmptyMatrix)
{
  TestEmptyVectorNonEmptyMatrix<Group1MatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group2MatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group3MatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group4MatrixAlias>();
}

TEST(VectorMatrix, NonEmptyVectorEmptyMatrix)
{
  TestNonEmptyVectorEmptyMatrix<Group1MatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group2MatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group3MatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group4MatrixAlias>();
}

TEST(VectorMatrix, EmptyVectorEmptyMatrix)
{
  TestEmptyVectorEmptyMatrix<Group1MatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group2MatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group3MatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<Group1MatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group2MatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group3MatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<Group1MatrixAlias>();
  TestMultipleVectorsSameSize<Group2MatrixAlias>();
  TestMultipleVectorsSameSize<Group3MatrixAlias>();
  TestMultipleVectorsSameSize<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleMatricesOneVector)
{
  TestMultipleMatricesOneVector<Group1MatrixAlias>();
  TestMultipleMatricesOneVector<Group2MatrixAlias>();
  TestMultipleMatricesOneVector<Group3MatrixAlias>();
  TestMultipleMatricesOneVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, MultipleMatricesDifferentRowsVector)
{
  TestMultipleMatricesDifferentRowsVector<Group1MatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group2MatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group3MatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, VectorSizeMatchesOneMatrixOnly)
{
  TestVectorSizeMatchesOneMatrixOnly<Group1MatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group2MatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group3MatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group4MatrixAlias>();
}

TEST(VectorMatrix, ConstVector)
{
  TestConstVector<Group1MatrixAlias>();
  TestConstVector<Group2MatrixAlias>();
  TestConstVector<Group3MatrixAlias>();
  TestConstVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, MutableVector)
{
  TestMutableVector<Group1MatrixAlias>();
  TestMutableVector<Group2MatrixAlias>();
  TestMutableVector<Group3MatrixAlias>();
  TestMutableVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, FunctionReusabilityWithVectors)
{
  TestFunctionReusabilityWithVectors<Group1MatrixAlias>();
  TestFunctionReusabilityWithVectors<Group2MatrixAlias>();
  TestFunctionReusabilityWithVectors<Group3MatrixAlias>();
  TestFunctionReusabilityWithVectors<Group4MatrixAlias>();
}

TEST(VectorMatrix, FunctionInvocationWithWrongSizedVector)
{
  TestFunctionInvocationWithWrongSizedVector<Group1MatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group2MatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group3MatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, ArraySupport)
{
  TestArraySupport<Group1MatrixAlias>();
  TestArraySupport<Group2MatrixAlias>();
  TestArraySupport<Group3MatrixAlias>();
  TestArraySupport<Group4MatrixAlias>();
}

TEST(VectorMatrix, MixedVectorColumnViewRowVariable)
{
  TestMixedVectorColumnViewRowVariable<Group1MatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group2MatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group3MatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group4MatrixAlias>();
}

TEST(VectorMatrix, IntegerVector)
{
  TestIntegerVector<Group1MatrixAlias>();
  TestIntegerVector<Group2MatrixAlias>();
  TestIntegerVector<Group3MatrixAlias>();
  TestIntegerVector<Group4MatrixAlias>();
}

TEST(VectorMatrix, FunctionWithConstSignature)
{
  TestFunctionWithConstSignature<Group1MatrixAlias>();
  TestFunctionWithConstSignature<Group2MatrixAlias>();
  TestFunctionWithConstSignature<Group3MatrixAlias>();
  TestFunctionWithConstSignature<Group4MatrixAlias>();
}
