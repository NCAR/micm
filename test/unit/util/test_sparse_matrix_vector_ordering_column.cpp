#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedColumnMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[12], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 16);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[16], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(), 5 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedColumnMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 8);
    EXPECT_EQ(matrix.AsVector()[8], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(SparseVectorCompressedColumnMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 8);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[8], 21);
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

TEST(SparseVectorCompressedColumnMatrix, Print)
{
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, PrintNonZero)
{
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

// Phase 2: BlockFunction infrastructure tests
TEST(SparseVectorCompressedColumnMatrix, ArrayFunction)
{
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, BlockViewReuse)
{
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionReusability)
{
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, TwoSparseMatricesDifferentStructure)
{
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseAndDenseMatrixFunction)
{
  // Only L=1 is compatible with standard-ordered dense Matrix
  testSparseAndDenseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>, micm::Matrix>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseAndVectorMatrixFunction)
{
  // Valid: Vector-ordered sparse with matching L vector matrix
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>, 1>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 2>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>, 3>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 4>();
  testSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<10>, 10>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleOrdering)
{
  // Invalid: L>1 vector-ordered sparse with L=1 standard-ordered dense
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, micm::Matrix>();
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>, micm::Matrix>();
  testIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, micm::Matrix>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleVectorOrdering)
{
  // Invalid: Vector-ordered sparse with mismatched L vector matrix
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 1>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 3>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 2>();
  testIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 10>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleSparseOrdering)
{
  // Invalid: Two sparse matrices with different L values
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
  testIncompatibleSparseOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, micm::SparseMatrixVectorOrderingCompressedSparseColumn<10>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedBlockDimensions)
{
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedElementDimensions)
{
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstSparseMatrixFunction)
{
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptySparseMatrixFunction)
{
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedBlocksAtInvocation)
{
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, WrongStructureAtInvocation)
{
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

// ============================================================================
// Vector Support Tests (using L=1 for simplicity)
// ============================================================================

TEST(SparseVectorCompressedColumnMatrix, VectorInSparseMatrixFunction)
{
  testVectorInSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorTooSmall)
{
  testVectorTooSmall<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorTooLarge)
{
  testVectorTooLarge<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptyVectorNonEmptySparseMatrix)
{
  testEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, NonEmptyVectorEmptySparseMatrix)
{
  testNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptyVectorEmptySparseMatrix)
{
  testEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleVectorsDifferentSizes)
{
  testMultipleVectorsDifferentSizes<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleVectorsSameSize)
{
  testMultipleVectorsSameSize<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesOneVector)
{
  testMultipleSparseMatricesOneVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  testMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  testVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstVectorSparse)
{
  testConstVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MutableVectorSparse)
{
  testMutableVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionReusabilityWithVectorsSparse)
{
  testFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  testFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, ArraySupportSparse)
{
  testArraySupportSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MixedVectorBlockViewBlockVariable)
{
  testMixedVectorBlockViewBlockVariable<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, IntegerVectorSparse)
{
  testIntegerVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}
