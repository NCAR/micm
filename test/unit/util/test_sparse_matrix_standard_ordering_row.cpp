#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

using StandardOrdering = micm::SparseMatrixStandardOrderingCompressedSparseRow;

TEST(SparseCompressedRowMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SingleBlockMatrix)
{
  {
    auto matrix = testSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 42;
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
}

TEST(SparseCompressedRowMatrix, ConstSingleBlockMatrix)
{
  {
    auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
}

TEST(SparseCompressedRowMatrix, MultiBlockMatrix)
{
  {
    auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(0, 2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 2, 1);
      EXPECT_EQ(elem, 12);
      matrix.AsVector()[elem] = 31;
      EXPECT_EQ(matrix.AsVector()[12], 31);
    }
  }
}

TEST(SparseCompressedRowMatrix, Print)
{
  testPrint<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, PrintNonZero)
{
  testPrintNonZero<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(3, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(2, 4); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(6, 7); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
}

TEST(SparseCompressedRowMatrix, ArrayFunction)
{
  testArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedBlockDimensions)
{
  testMismatchedBlockDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedElementDimensions)
{
  testMismatchedElementDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, BlockViewReuse)
{
  testBlockViewReuse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, FunctionReusability)
{
  testFunctionReusability<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, TwoSparseMatricesDifferentStructure)
{
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SparseAndDenseMatrixFunction)
{
  testSparseAndDenseMatrixFunction<micm::SparseMatrix, StandardOrdering, micm::Matrix>();
}

TEST(SparseCompressedRowMatrix, ConstSparseMatrixFunction)
{
  testConstSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptySparseMatrixFunction)
{
  testEmptySparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedBlocksAtInvocation)
{
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, WrongStructureAtInvocation)
{
  testWrongStructureAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

// ============================================================================
// Vector Support Tests
// ============================================================================

TEST(SparseCompressedRowMatrix, VectorInSparseMatrixFunction)
{
  testVectorInSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorTooSmall)
{
  testVectorTooSmall<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorTooLarge)
{
  testVectorTooLarge<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptyVectorNonEmptySparseMatrix)
{
  testEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, NonEmptyVectorEmptySparseMatrix)
{
  testNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptyVectorEmptySparseMatrix)
{
  testEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleVectorsDifferentSizes)
{
  testMultipleVectorsDifferentSizes<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleVectorsSameSize)
{
  testMultipleVectorsSameSize<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesOneVector)
{
  testMultipleSparseMatricesOneVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  testMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  testVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ConstVectorSparse)
{
  testConstVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MutableVectorSparse)
{
  testMutableVectorSparse<micm::SparseMatrix, StandardOrdering>();
}
TEST(SparseCompressedRowMatrix, FunctionWithConstSignatureSparse)
{
  testFunctionWithConstSignatureSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, GetBlockViewByVectorIndex)
{
  testGetBlockViewByVectorIndex<micm::SparseMatrix, StandardOrdering>();
}
TEST(SparseCompressedRowMatrix, FunctionReusabilityWithVectorsSparse)
{
  testFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  testFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ArraySupportSparse)
{
  testArraySupportSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MixedVectorBlockViewBlockVariable)
{
  testMixedVectorBlockViewBlockVariable<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, IntegerVectorSparse)
{
  testIntegerVectorSparse<micm::SparseMatrix, StandardOrdering>();
}