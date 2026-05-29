#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

TEST(SparseVectorCompressedColumnMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, SetScalar)
{
  TestSetScalar<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestAddToDiagonal<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SingleBlockMatrix)
{
  auto matrix = TestSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();

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
  auto matrix = TestConstSingleBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

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
  auto matrix = TestMultiBlockMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();

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
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestPrint<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestPrintNonZero<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

// Phase 2: BlockFunction infrastructure tests
TEST(SparseVectorCompressedColumnMatrix, ArrayFunction)
{
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMultiMatrixArrayFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMultipleTemporaries<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, BlockViewReuse)
{
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestBlockViewReuse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionReusability)
{
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestFunctionReusability<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, TwoSparseMatricesDifferentStructure)
{
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseAndDenseMatrixFunction)
{
  // Only L=1 is compatible with standard-ordered dense Matrix
  TestSparseAndDenseMatrixFunction<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>,
      micm::Matrix>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseAndVectorMatrixFunction)
{
  // Valid: Vector-ordered sparse with matching L vector matrix
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>, 1>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 2>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>, 3>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 4>();
  TestSparseAndVectorMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<10>, 10>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleOrdering)
{
  // Invalid: L>1 vector-ordered sparse with L=1 standard-ordered dense
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, micm::Matrix>();
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>, micm::Matrix>();
  TestIncompatibleOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, micm::Matrix>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleVectorOrdering)
{
  // Invalid: Vector-ordered sparse with mismatched L vector matrix
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 1>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>, 3>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 2>();
  TestIncompatibleVectorOrdering<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>, 10>();
}

TEST(SparseVectorCompressedColumnMatrix, IncompatibleSparseOrdering)
{
  // Invalid: Two sparse matrices with different L values
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
  TestIncompatibleSparseOrdering<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<10>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedBlockDimensions)
{
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMismatchedBlockDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedElementDimensions)
{
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMismatchedElementDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestWrongMatrixDimensions<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstSparseMatrixFunction)
{
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestConstSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptySparseMatrixFunction)
{
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestEmptySparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MismatchedBlocksAtInvocation)
{
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, WrongStructureAtInvocation)
{
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestWrongStructureAtInvocation<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

// ============================================================================
// Vector Support Tests (using L=1 for simplicity)
// ============================================================================

TEST(SparseVectorCompressedColumnMatrix, VectorInSparseMatrixFunction)
{
  TestVectorInSparseMatrixFunction<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptyVectorNonEmptySparseMatrix)
{
  TestEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, NonEmptyVectorEmptySparseMatrix)
{
  TestNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, EmptyVectorEmptySparseMatrix)
{
  TestEmptyVectorEmptySparseMatrix<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesOneVector)
{
  TestMultipleSparseMatricesOneVector<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  TestMultipleSparseMatricesDifferentBlocksVector<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  TestVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstVectorSparse)
{
  TestConstVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MutableVectorSparse)
{
  TestMutableVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionReusabilityWithVectorsSparse)
{
  TestFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  TestFunctionInvocationWithWrongSizedVectorSparse<
      micm::SparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, ArraySupportSparse)
{
  TestArraySupportSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, MixedVectorBlockViewBlockVariable)
{
  TestMixedVectorBlockViewBlockVariable<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, IntegerVectorSparse)
{
  TestIntegerVectorSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
}

TEST(SparseVectorCompressedColumnMatrix, FunctionWithConstSignatureSparse)
{
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, GetBlockViewByVectorIndex)
{
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}
