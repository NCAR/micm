#include "test_sparse_matrix_policy.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

using StandardOrdering = micm::SparseMatrixStandardOrderingCompressedSparseColumn;

TEST(SparseCompressedColumnMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SetScalar)
{
  TestSetScalar<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SingleBlockMatrix)
{
  {
    auto matrix = TestSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();

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
    auto matrix = TestConstSingleBlockMatrix<micm::SparseMatrix, StandardOrdering>();
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
    auto matrix = TestMultiBlockMatrix<micm::SparseMatrix, StandardOrdering>();

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
  auto matrix = TestPrint<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseMatrixBuilder, BadConfiguration)
{
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(3, 0); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(2, 4); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
  EXPECT_THROW(
      try { auto builder = micm::SparseMatrix<double>::Create(3).WithElement(6, 7); } catch (micm::MicmException& e) {
        EXPECT_EQ(e.code_, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE);
        throw;
      },
      micm::MicmException);
}

TEST(SparseCompressedColumnMatrix, ArrayFunction)
{
  TestArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedBlockDimensions)
{
  TestMismatchedBlockDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedElementDimensions)
{
  TestMismatchedElementDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, WrongMatrixDimensions)
{
  TestWrongMatrixDimensions<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, BlockViewReuse)
{
  TestBlockViewReuse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionReusability)
{
  TestFunctionReusability<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, TwoSparseMatricesDifferentStructure)
{
  TestTwoSparseMatricesDifferentStructure<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SparseAndDenseMatrixFunction)
{
  TestSparseAndDenseMatrixFunction<micm::SparseMatrix, StandardOrdering, micm::Matrix>();
}

TEST(SparseCompressedColumnMatrix, ConstSparseMatrixFunction)
{
  TestConstSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptySparseMatrixFunction)
{
  TestEmptySparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MismatchedBlocksAtInvocation)
{
  TestMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, WrongStructureAtInvocation)
{
  TestWrongStructureAtInvocation<micm::SparseMatrix, StandardOrdering>();
}

// ============================================================================
// Vector Support Tests
// ============================================================================

TEST(SparseCompressedColumnMatrix, VectorInSparseMatrixFunction)
{
  TestVectorInSparseMatrixFunction<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptyVectorNonEmptySparseMatrix)
{
  TestEmptyVectorNonEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, NonEmptyVectorEmptySparseMatrix)
{
  TestNonEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, EmptyVectorEmptySparseMatrix)
{
  TestEmptyVectorEmptySparseMatrix<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesOneVector)
{
  TestMultipleSparseMatricesOneVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  TestMultipleSparseMatricesDifferentBlocksVector<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  TestVectorSizeMatchesOneSparseMatrixOnly<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ConstVectorSparse)
{
  TestConstVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MutableVectorSparse)
{
  TestMutableVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionReusabilityWithVectorsSparse)
{
  TestFunctionReusabilityWithVectorsSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  TestFunctionInvocationWithWrongSizedVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, ArraySupportSparse)
{
  TestArraySupportSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, MixedVectorBlockViewBlockVariable)
{
  TestMixedVectorBlockViewBlockVariable<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, IntegerVectorSparse)
{
  TestIntegerVectorSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, FunctionWithConstSignatureSparse)
{
  TestFunctionWithConstSignatureSparse<micm::SparseMatrix, StandardOrdering>();
}

TEST(SparseCompressedColumnMatrix, GetBlockViewByVectorIndex)
{
  TestGetBlockViewByVectorIndex<micm::SparseMatrix, StandardOrdering>();
}
