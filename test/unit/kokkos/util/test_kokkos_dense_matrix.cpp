#include "../../util/test_matrix_policy.hpp"

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

template<class T>
using Group1KokkosMatrixAlias = micm::KokkosDenseMatrix<T, 1>;
template<class T>
using Group2KokkosMatrixAlias = micm::KokkosDenseMatrix<T, 2>;
template<class T>
using Group3KokkosMatrixAlias = micm::KokkosDenseMatrix<T, 3>;
template<class T>
using Group4KokkosMatrixAlias = micm::KokkosDenseMatrix<T, 4>;

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
  constexpr micm::Index kGroups = (3 + MICM_DEFAULT_VECTOR_SIZE - 1) / MICM_DEFAULT_VECTOR_SIZE;
  EXPECT_EQ(matrix.AsVector().size(), kGroups * MICM_DEFAULT_VECTOR_SIZE * 4);
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
  for (auto& elem : matrix.AsVector())
  {
    elem = 0.0;
  }
  EXPECT_EQ(matrix[0][0], 0.0);

  matrix.CopyToHost();
  EXPECT_EQ(matrix[0][0], 1.0);
  EXPECT_EQ(matrix[0][1], 2.0);
  EXPECT_EQ(matrix[1][0], 3.0);
  EXPECT_EQ(matrix[1][1], 4.0);
}

TEST(KokkosDenseMatrix, SmallVectorMatrix)
{
  auto matrix = TestSmallMatrix<Group2KokkosMatrixAlias>();

  std::vector<micm::Real>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 2);
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(data[0], static_cast<micm::Real>(41.2));
  EXPECT_EQ(data[2 * 5 + 0 + 2 * 4], static_cast<micm::Real>(102.3));
  EXPECT_EQ(data[1 + 2 * 3], static_cast<micm::Real>(64.7));
}

TEST(KokkosDenseMatrix, SmallConstVectorMatrix)
{
  auto matrix = TestSmallConstMatrix<Group4KokkosMatrixAlias>();

  const std::vector<micm::Real>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 4 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 1);
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(data[0], static_cast<micm::Real>(41.2));
  EXPECT_EQ(data[2 + 4 * 4], static_cast<micm::Real>(102.3));
  EXPECT_EQ(data[1 + 4 * 3], static_cast<micm::Real>(64.7));
}

TEST(KokkosDenseMatrix, InitializeVectorMatrix)
{
  TestInializeMatrix<Group1KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, InitializeConstVectorMatrix)
{
  TestInializeConstMatrix<Group2KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, LoopOverVectorMatrix)
{
  TestLoopOverMatrix<Group2KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, LoopOverConstVectorMatrix)
{
  TestLoopOverConstMatrix<Group1KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, Strides)
{
  auto matrix3vec = TestStrides<Group3KokkosMatrixAlias>();
  EXPECT_EQ(matrix3vec.RowStride(), 1);
  EXPECT_EQ(matrix3vec.ColumnStride(), 3);
  auto matrix4vec = TestStrides<Group4KokkosMatrixAlias>();
  EXPECT_EQ(matrix4vec.RowStride(), 1);
  EXPECT_EQ(matrix4vec.ColumnStride(), 4);
}

TEST(KokkosDenseMatrix, ConversionToVector)
{
  TestConversionToVector<Group3KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ConstConversionToVector)
{
  TestConstConversionToVector<Group1KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ConversionFromVector)
{
  TestConversionFromVector<Group2KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, AssignmentFromVector)
{
  TestAssignmentFromVector<Group2KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, Axpy)
{
  TestAxpy<Group1KokkosMatrixAlias>();
  TestAxpy<Group2KokkosMatrixAlias>();
  TestAxpy<Group3KokkosMatrixAlias>();
  TestAxpy<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, SetScaler)
{
  TestSetScalar<Group1KokkosMatrixAlias>();
  TestSetScalar<Group2KokkosMatrixAlias>();
  TestSetScalar<Group3KokkosMatrixAlias>();
  TestSetScalar<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, Max)
{
  TestMax<Group1KokkosMatrixAlias>();
  TestMax<Group2KokkosMatrixAlias>();
  TestMax<Group3KokkosMatrixAlias>();
  TestMax<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, Min)
{
  TestMin<Group1KokkosMatrixAlias>();
  TestMin<Group2KokkosMatrixAlias>();
  TestMin<Group3KokkosMatrixAlias>();
  TestMin<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, Print)
{
  TestPrint<Group1KokkosMatrixAlias>();
  TestPrint<Group2KokkosMatrixAlias>();
  TestPrint<Group3KokkosMatrixAlias>();
  TestPrint<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ForEach)
{
  TestForEach<Group1KokkosMatrixAlias>();
  TestForEach<Group2KokkosMatrixAlias>();
  TestForEach<Group3KokkosMatrixAlias>();
  TestForEach<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ArrayFunction)
{
  TestArrayFunction<Group1KokkosMatrixAlias>();
  TestArrayFunction<Group2KokkosMatrixAlias>();
  TestArrayFunction<Group3KokkosMatrixAlias>();
  TestArrayFunction<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<Group1KokkosMatrixAlias>();
  TestMultiMatrixArrayFunction<Group2KokkosMatrixAlias>();
  TestMultiMatrixArrayFunction<Group3KokkosMatrixAlias>();
  TestMultiMatrixArrayFunction<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MismatchedRowDimensions)
{
  TestMismatchedRowDimensions<Group1KokkosMatrixAlias>();
  TestMismatchedRowDimensions<Group2KokkosMatrixAlias>();
  TestMismatchedRowDimensions<Group3KokkosMatrixAlias>();
  TestMismatchedRowDimensions<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MismatchedColumnDimensions)
{
  TestMismatchedColumnDimensions<Group1KokkosMatrixAlias>();
  TestMismatchedColumnDimensions<Group2KokkosMatrixAlias>();
  TestMismatchedColumnDimensions<Group3KokkosMatrixAlias>();
  TestMismatchedColumnDimensions<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<Group1KokkosMatrixAlias>();
  TestWrongMatrixDimensions<Group2KokkosMatrixAlias>();
  TestWrongMatrixDimensions<Group3KokkosMatrixAlias>();
  TestWrongMatrixDimensions<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<Group1KokkosMatrixAlias>();
  TestMultipleTemporaries<Group2KokkosMatrixAlias>();
  TestMultipleTemporaries<Group3KokkosMatrixAlias>();
  TestMultipleTemporaries<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ColumnViewReuse)
{
  TestColumnViewReuse<Group1KokkosMatrixAlias>();
  TestColumnViewReuse<Group2KokkosMatrixAlias>();
  TestColumnViewReuse<Group3KokkosMatrixAlias>();
  TestColumnViewReuse<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, FunctionReusability)
{
  TestFunctionReusability<Group1KokkosMatrixAlias>();
  TestFunctionReusability<Group2KokkosMatrixAlias>();
  TestFunctionReusability<Group3KokkosMatrixAlias>();
  TestFunctionReusability<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ConstMatrixFunction)
{
  TestConstMatrixFunction<Group1KokkosMatrixAlias>();
  TestConstMatrixFunction<Group2KokkosMatrixAlias>();
  TestConstMatrixFunction<Group3KokkosMatrixAlias>();
  TestConstMatrixFunction<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, EmptyMatrixFunction)
{
  TestEmptyMatrixFunction<Group1KokkosMatrixAlias>();
  TestEmptyMatrixFunction<Group2KokkosMatrixAlias>();
  TestEmptyMatrixFunction<Group3KokkosMatrixAlias>();
  TestEmptyMatrixFunction<Group4KokkosMatrixAlias>();
}

// Flexible row count Tests
TEST(KokkosDenseMatrix, MultiMatrixDifferentRowsFromCreation)
{
  TestMultiMatrixDifferentRowsFromCreation<Group1KokkosMatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group2KokkosMatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group3KokkosMatrixAlias>();
  TestMultiMatrixDifferentRowsFromCreation<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MatrixVectorDifferentRowsFromCreation)
{
  TestMatrixVectorDifferentRowsFromCreation<Group1KokkosMatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group2KokkosMatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group3KokkosMatrixAlias>();
  TestMatrixVectorDifferentRowsFromCreation<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MismatchedRowsAtInvocation)
{
  TestMismatchedRowsAtInvocation<Group1KokkosMatrixAlias>();
  TestMismatchedRowsAtInvocation<Group2KokkosMatrixAlias>();
  TestMismatchedRowsAtInvocation<Group3KokkosMatrixAlias>();
  TestMismatchedRowsAtInvocation<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleMatricesMismatchedRowsAtInvocation)
{
  TestMultipleMatricesMismatchedRowsAtInvocation<Group1KokkosMatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group2KokkosMatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group3KokkosMatrixAlias>();
  TestMultipleMatricesMismatchedRowsAtInvocation<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, WrongColumnCountAtInvocation)
{
  TestWrongColumnCountAtInvocation<Group1KokkosMatrixAlias>();
  TestWrongColumnCountAtInvocation<Group2KokkosMatrixAlias>();
  TestWrongColumnCountAtInvocation<Group3KokkosMatrixAlias>();
  TestWrongColumnCountAtInvocation<Group4KokkosMatrixAlias>();
}

// Vector support Tests
TEST(KokkosDenseMatrix, VectorInMatrixFunction)
{
  TestVectorInMatrixFunction<Group1KokkosMatrixAlias>();
  TestVectorInMatrixFunction<Group2KokkosMatrixAlias>();
  TestVectorInMatrixFunction<Group3KokkosMatrixAlias>();
  TestVectorInMatrixFunction<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, VectorTooSmall)
{
  TestVectorTooSmall<Group1KokkosMatrixAlias>();
  TestVectorTooSmall<Group2KokkosMatrixAlias>();
  TestVectorTooSmall<Group3KokkosMatrixAlias>();
  TestVectorTooSmall<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, VectorTooLarge)
{
  TestVectorTooLarge<Group1KokkosMatrixAlias>();
  TestVectorTooLarge<Group2KokkosMatrixAlias>();
  TestVectorTooLarge<Group3KokkosMatrixAlias>();
  TestVectorTooLarge<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, EmptyVectorNonEmptyMatrix)
{
  TestEmptyVectorNonEmptyMatrix<Group1KokkosMatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group2KokkosMatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group3KokkosMatrixAlias>();
  TestEmptyVectorNonEmptyMatrix<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, NonEmptyVectorEmptyMatrix)
{
  TestNonEmptyVectorEmptyMatrix<Group1KokkosMatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group2KokkosMatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group3KokkosMatrixAlias>();
  TestNonEmptyVectorEmptyMatrix<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, EmptyVectorEmptyMatrix)
{
  TestEmptyVectorEmptyMatrix<Group1KokkosMatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group2KokkosMatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group3KokkosMatrixAlias>();
  TestEmptyVectorEmptyMatrix<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<Group1KokkosMatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group2KokkosMatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group3KokkosMatrixAlias>();
  TestMultipleVectorsDifferentSizes<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<Group1KokkosMatrixAlias>();
  TestMultipleVectorsSameSize<Group2KokkosMatrixAlias>();
  TestMultipleVectorsSameSize<Group3KokkosMatrixAlias>();
  TestMultipleVectorsSameSize<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleMatricesOneVector)
{
  TestMultipleMatricesOneVector<Group1KokkosMatrixAlias>();
  TestMultipleMatricesOneVector<Group2KokkosMatrixAlias>();
  TestMultipleMatricesOneVector<Group3KokkosMatrixAlias>();
  TestMultipleMatricesOneVector<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MultipleMatricesDifferentRowsVector)
{
  TestMultipleMatricesDifferentRowsVector<Group1KokkosMatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group2KokkosMatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group3KokkosMatrixAlias>();
  TestMultipleMatricesDifferentRowsVector<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, VectorSizeMatchesOneMatrixOnly)
{
  TestVectorSizeMatchesOneMatrixOnly<Group1KokkosMatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group2KokkosMatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group3KokkosMatrixAlias>();
  TestVectorSizeMatchesOneMatrixOnly<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ConstVector)
{
  TestConstVector<Group1KokkosMatrixAlias>();
  TestConstVector<Group2KokkosMatrixAlias>();
  TestConstVector<Group3KokkosMatrixAlias>();
  TestConstVector<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MutableVector)
{
  TestMutableVector<Group1KokkosMatrixAlias>();
  TestMutableVector<Group2KokkosMatrixAlias>();
  TestMutableVector<Group3KokkosMatrixAlias>();
  TestMutableVector<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, FunctionReusabilityWithVectors)
{
  TestFunctionReusabilityWithVectors<Group1KokkosMatrixAlias>();
  TestFunctionReusabilityWithVectors<Group2KokkosMatrixAlias>();
  TestFunctionReusabilityWithVectors<Group3KokkosMatrixAlias>();
  TestFunctionReusabilityWithVectors<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, FunctionInvocationWithWrongSizedVector)
{
  TestFunctionInvocationWithWrongSizedVector<Group1KokkosMatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group2KokkosMatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group3KokkosMatrixAlias>();
  TestFunctionInvocationWithWrongSizedVector<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, MixedVectorColumnViewRowVariable)
{
  TestMixedVectorColumnViewRowVariable<Group1KokkosMatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group2KokkosMatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group3KokkosMatrixAlias>();
  TestMixedVectorColumnViewRowVariable<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, IntegerVector)
{
  TestIntegerVector<Group1KokkosMatrixAlias>();
  TestIntegerVector<Group2KokkosMatrixAlias>();
  TestIntegerVector<Group3KokkosMatrixAlias>();
  TestIntegerVector<Group4KokkosMatrixAlias>();
}

#ifndef KOKKOS_ENABLE_CUDA
// This test attemps to wrap the function in std::function which isn't allowed
// with CUDA
TEST(KokkosDenseMatrix, FunctionWithConstSignature)
{
  TestFunctionWithConstSignature<Group1KokkosMatrixAlias>();
  TestFunctionWithConstSignature<Group2KokkosMatrixAlias>();
  TestFunctionWithConstSignature<Group3KokkosMatrixAlias>();
  TestFunctionWithConstSignature<Group4KokkosMatrixAlias>();
}
#endif

TEST(KokkosDenseMatrix, TestFill)
{
  TestFill<Group1KokkosMatrixAlias>();
  TestFill<Group2KokkosMatrixAlias>();
  TestFill<Group3KokkosMatrixAlias>();
  TestFill<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, TestCopy)
{
  TestCopy<Group1KokkosMatrixAlias>();
  TestCopy<Group2KokkosMatrixAlias>();
  TestCopy<Group3KokkosMatrixAlias>();
  TestCopy<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ReduceSum)
{
  TestReduceSum<Group1KokkosMatrixAlias>();
  TestReduceSum<Group2KokkosMatrixAlias>();
  TestReduceSum<Group3KokkosMatrixAlias>();
  TestReduceSum<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ReduceMax)
{
  TestReduceMax<Group1KokkosMatrixAlias>();
  TestReduceMax<Group2KokkosMatrixAlias>();
  TestReduceMax<Group3KokkosMatrixAlias>();
  TestReduceMax<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ReduceLOr)
{
  TestReduceLOr<Group1KokkosMatrixAlias>();
  TestReduceLOr<Group2KokkosMatrixAlias>();
  TestReduceLOr<Group3KokkosMatrixAlias>();
  TestReduceLOr<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ReduceLAnd)
{
  TestReduceLAnd<Group1KokkosMatrixAlias>();
  TestReduceLAnd<Group2KokkosMatrixAlias>();
  TestReduceLAnd<Group3KokkosMatrixAlias>();
  TestReduceLAnd<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, ReduceStrict)
{
  TestReduceStrict<Group1KokkosMatrixAlias>();
  TestReduceStrict<Group2KokkosMatrixAlias>();
  TestReduceStrict<Group3KokkosMatrixAlias>();
  TestReduceStrict<Group4KokkosMatrixAlias>();
}

TEST(KokkosDenseMatrix, VectorCapture)
{
  TestVectorCapture<Group1KokkosMatrixAlias>();
  TestVectorCapture<Group2KokkosMatrixAlias>();
  TestVectorCapture<Group3KokkosMatrixAlias>();
  TestVectorCapture<Group4KokkosMatrixAlias>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
