#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedRowMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, SetScalar)
{
  TestSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SingleBlockMatrix)
{
  auto matrix = TestSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();

  {
    micm::Index elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[16], 42);
  }
  {
    micm::Index elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[12], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, ConstSingleBlockMatrix)
{
  auto matrix = TestConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    micm::Index elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 42);
  }
  {
    micm::Index elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, MultiBlockMatrix)
{
  auto matrix = TestMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    micm::Index elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 6);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  {
    micm::Index elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 14);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[14], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}

TEST(SparseVectorCompressedRowMatrix, Print)
{
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

// Phase 2: BlockFunction infrastructure Tests
TEST(SparseVectorCompressedRowMatrix, ArrayFunction)
{
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, BlockViewReuse)
{
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionReusability)
{
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, TwoSparseMatricesDifferentStructure)
{
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SparseAndDenseMatrixFunction)
{
  // Only L=1 is compatible with standard-ordered dense Matrix
  TestSparseAndDenseMatrixFunction<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>,
      micm::Matrix>();
}

TEST(SparseVectorCompressedRowMatrix, SparseAndVectorMatrixFunction)
{
  // Valid: Vector-ordered sparse with matching L vector matrix
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>, 1>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 2>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>, 3>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 4>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<10>, 10>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleOrdering)
{
  // Invalid: L>1 vector-ordered sparse with L=1 standard-ordered dense
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, micm::Matrix>();
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>, micm::Matrix>();
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, micm::Matrix>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleVectorOrdering)
{
  // Invalid: Vector-ordered sparse with mismatched L vector matrix
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 1>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 3>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 2>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 10>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleSparseOrdering)
{
  // Invalid: Two sparse matrices with different L values
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<2>,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<3>,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<4>,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<10>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedBlockDimensions)
{
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedElementDimensions)
{
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstSparseMatrixFunction)
{
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptySparseMatrixFunction)
{
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedBlocksAtInvocation)
{
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, WrongStructureAtInvocation)
{
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

// ============================================================================
// Vector Support Tests (using L=1 for simplicity)
// ============================================================================

TEST(SparseVectorCompressedRowMatrix, VectorInSparseMatrixFunction)
{
  TestVectorInSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptyVectorNonEmptySparseMatrix)
{
  TestEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, NonEmptyVectorEmptySparseMatrix)
{
  TestNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptyVectorEmptySparseMatrix)
{
  TestEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesOneVector)
{
  TestMultipleSparseMatricesOneVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  TestMultipleSparseMatricesDifferentBlocksVector<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  TestVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstVectorSparse)
{
  TestConstVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MutableVectorSparse)
{
  TestMutableVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionReusabilityWithVectorsSparse)
{
  TestFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  TestFunctionInvocationWithWrongSizedVectorSparse<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, ArraySupportSparse)
{
  TestArraySupportSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MixedVectorBlockViewBlockVariable)
{
  TestMixedVectorBlockViewBlockVariable<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, IntegerVectorSparse)
{
  TestIntegerVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionWithConstSignatureSparse)
{
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, GetBlockViewByVectorIndex)
{
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}
