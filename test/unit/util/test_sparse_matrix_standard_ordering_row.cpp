#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using StandardOrdering = micm::SparseMatrixStandardOrderingCompressedSparseRow;

TEST(SparseCompressedRowMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SetScalar)
{
  TestSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SingleBlockMatrix)
{
  {
    auto matrix = TestSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      micm::Index elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      matrix.AsVector()[elem] = 42;
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      micm::Index elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
}

TEST(SparseCompressedRowMatrix, ConstSingleBlockMatrix)
{
  {
    auto matrix = TestConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
    {
      micm::Index elem = matrix.VectorIndex(3, 2);
      EXPECT_EQ(elem, 4);
      EXPECT_EQ(matrix.AsVector()[4], 42);
    }
    {
      micm::Index elem = matrix.VectorIndex(2, 3);
      EXPECT_EQ(elem, 3);
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
  }
}

TEST(SparseCompressedRowMatrix, MultiBlockMatrix)
{
  {
    auto matrix = TestMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

    {
      micm::Index elem = matrix.VectorIndex(0, 2, 3);
      EXPECT_EQ(elem, 3);
      matrix.AsVector()[elem] = 21;
      EXPECT_EQ(matrix.AsVector()[3], 21);
    }
    {
      micm::Index elem = matrix.VectorIndex(2, 2, 1);
      EXPECT_EQ(elem, 12);
      matrix.AsVector()[elem] = 31;
      EXPECT_EQ(matrix.AsVector()[12], 31);
    }
  }
}

TEST(SparseCompressedRowMatrix, Print)
{
  TestPrint<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<micm::Real>::Create(3).WithElement(3, 0); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<micm::Real>::Create(3).WithElement(2, 4); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<micm::Real>::Create(3).WithElement(6, 7); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
}

TEST(SparseCompressedRowMatrix, ArrayFunction)
{
  TestArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedBlockDimensions)
{
  TestMismatchedBlockDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedElementDimensions)
{
  TestMismatchedElementDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, BlockViewReuse)
{
  TestBlockViewReuse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, FunctionReusability)
{
  TestFunctionReusability<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, TwoSparseMatricesDifferentStructure)
{
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SparseAndDenseMatrixFunction)
{
  TestSparseAndDenseMatrixFunction<micm::SparseMatrix, StandardOrdering, micm::Matrix>();
}

TEST(SparseCompressedRowMatrix, ConstSparseMatrixFunction)
{
  TestConstSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptySparseMatrixFunction)
{
  TestEmptySparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MismatchedBlocksAtInvocation)
{
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, WrongStructureAtInvocation)
{
  TestWrongStructureAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

// ============================================================================
// Vector Support Tests
// ============================================================================

TEST(SparseCompressedRowMatrix, VectorInSparseMatrixFunction)
{
  TestVectorInSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptyVectorNonEmptySparseMatrix)
{
  TestEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, NonEmptyVectorEmptySparseMatrix)
{
  TestNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, EmptyVectorEmptySparseMatrix)
{
  TestEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesOneVector)
{
  TestMultipleSparseMatricesOneVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  TestMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  TestVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ConstVectorSparse)
{
  TestConstVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MutableVectorSparse)
{
  TestMutableVectorSparse<micm::SparseMatrix, StandardOrdering>();
}
TEST(SparseCompressedRowMatrix, FunctionWithConstSignatureSparse)
{
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, GetBlockViewByVectorIndex)
{
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, StandardOrdering>();
}
TEST(SparseCompressedRowMatrix, FunctionReusabilityWithVectorsSparse)
{
  TestFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  TestFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, ArraySupportSparse)
{
  TestArraySupportSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, MixedVectorBlockViewBlockVariable)
{
  TestMixedVectorBlockViewBlockVariable<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedRowMatrix, IntegerVectorSparse)
{
  TestIntegerVectorSparse<micm::SparseMatrix, StandardOrdering>();
}