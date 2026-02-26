#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

using StandardOrdering = micm::SparseMatrixStandardOrderingCompressedSparseColumn;

TEST(SparseCompressedColumnMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SetScalar)
{
  testSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SingleBlockMatrix)
{
  {
    auto matrix = testSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 42;
      EXPECT_EQ(matrix.AsVector()[3], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
  }
}

TEST(SparseCompressedColumnMatrix, ConstSingleBlockMatrix)
{
  {
    auto matrix = testConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
    {
      std::size_t elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 3);
      EXPECT_EQ(matrix.AsVector()[3], 42);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 4);
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
  }
}

TEST(SparseCompressedColumnMatrix, MultiBlockMatrix)
{
  {
    auto matrix = testMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      std::size_t elem = matrix.VectorIndex(0, 2, 3);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[4], 21);
    }
    {
      std::size_t elem = matrix.VectorIndex(2, 2, 1);
      EXPECT_EQ(elem, 12);
      matrix.AsVector()[elem] = 31;
      EXPECT_EQ(matrix.AsVector()[12], 31);
    }
  }
}

TEST(SparseCompressedColumnMatrix, Print)
{
  auto matrix = testPrint<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, PrintNonZero)
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

TEST(SparseCompressedColumnMatrix, ArrayFunction)
{
  testArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedBlockDimensions)
{
  testMismatchedBlockDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedElementDimensions)
{
  testMismatchedElementDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, BlockViewReuse)
{
  testBlockViewReuse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionReusability)
{
  testFunctionReusability<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, TwoSparseMatricesDifferentStructure)
{
  testTwoSparseMatricesDifferentStructure<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SparseAndDenseMatrixFunction)
{
  testSparseAndDenseMatrixFunction<micm::SparseMatrix, StandardOrdering, micm::Matrix>();
}

TEST(SparseCompressedColumnMatrix, ConstSparseMatrixFunction)
{
  testConstSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptySparseMatrixFunction)
{
  testEmptySparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedBlocksAtInvocation)
{
  testMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, WrongStructureAtInvocation)
{
  testWrongStructureAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

// ============================================================================
// Vector Support Tests
// ============================================================================

TEST(SparseCompressedColumnMatrix, VectorInSparseMatrixFunction)
{
  testVectorInSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorTooSmall)
{
  testVectorTooSmall<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorTooLarge)
{
  testVectorTooLarge<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptyVectorNonEmptySparseMatrix)
{
  testEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, NonEmptyVectorEmptySparseMatrix)
{
  testNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptyVectorEmptySparseMatrix)
{
  testEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleVectorsDifferentSizes)
{
  testMultipleVectorsDifferentSizes<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleVectorsSameSize)
{
  testMultipleVectorsSameSize<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesOneVector)
{
  testMultipleSparseMatricesOneVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  testMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  testVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ConstVectorSparse)
{
  testConstVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MutableVectorSparse)
{
  testMutableVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionReusabilityWithVectorsSparse)
{
  testFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  testFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ArraySupportSparse)
{
  testArraySupportSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MixedVectorBlockViewBlockVariable)
{
  testMixedVectorBlockViewBlockVariable<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, IntegerVectorSparse)
{
  testIntegerVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionWithConstSignatureSparse)
{
  testFunctionWithConstSignatureSparse<micm::SparseMatrix, StandardOrdering>();
}
