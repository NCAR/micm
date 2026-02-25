#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedRowMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
}

TEST(SparseVectorCompressedRowMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[16], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
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
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedRowMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 6);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[6], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
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
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, PrintNonZero)
{
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

// Phase 2: BlockFunction infrastructure tests
TEST(SparseVectorCompressedRowMatrix, ArrayFunction)
{
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, BlockViewReuse)
{
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionReusability)
{
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, TwoSparseMatricesDifferentStructure)
{
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SparseAndDenseMatrixFunction)
{
  // Only L=1 is compatible with standard-ordered dense Matrix
  testSparseAndDenseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>, micm::Matrix>();
}

TEST(SparseVectorCompressedRowMatrix, SparseAndVectorMatrixFunction)
{
  // Valid: Vector-ordered sparse with matching L vector matrix
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>, 1>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 2>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>, 3>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 4>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<10>, 10>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleOrdering)
{
  // Invalid: L>1 vector-ordered sparse with L=1 standard-ordered dense
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, micm::Matrix>();
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>, micm::Matrix>();
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, micm::Matrix>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleVectorOrdering)
{
  // Invalid: Vector-ordered sparse with mismatched L vector matrix
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 1>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, 3>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 2>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, 10>();
}

TEST(SparseVectorCompressedRowMatrix, IncompatibleSparseOrdering)
{
  // Invalid: Two sparse matrices with different L values
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>, micm::SparseMatrixVectorOrderingCompressedSparseRow<10>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedBlockDimensions)
{
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedElementDimensions)
{
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstSparseMatrixFunction)
{
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptySparseMatrixFunction)
{
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MismatchedBlocksAtInvocation)
{
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

TEST(SparseVectorCompressedRowMatrix, WrongStructureAtInvocation)
{
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<3>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>();
}

// ============================================================================
// Vector Support Tests (using L=1 for simplicity)
// ============================================================================

TEST(SparseVectorCompressedRowMatrix, VectorInSparseMatrixFunction)
{
  testVectorInSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorTooSmall)
{
  testVectorTooSmall<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorTooLarge)
{
  testVectorTooLarge<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptyVectorNonEmptySparseMatrix)
{
  testEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, NonEmptyVectorEmptySparseMatrix)
{
  testNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, EmptyVectorEmptySparseMatrix)
{
  testEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleVectorsDifferentSizes)
{
  testMultipleVectorsDifferentSizes<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleVectorsSameSize)
{
  testMultipleVectorsSameSize<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesOneVector)
{
  testMultipleSparseMatricesOneVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  testMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  testVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, ConstVectorSparse)
{
  testConstVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MutableVectorSparse)
{
  testMutableVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionReusabilityWithVectorsSparse)
{
  testFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  testFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, ArraySupportSparse)
{
  testArraySupportSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, MixedVectorBlockViewBlockVariable)
{
  testMixedVectorBlockViewBlockVariable<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}

TEST(SparseVectorCompressedRowMatrix, IntegerVectorSparse)
{
  testIntegerVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>();
}
