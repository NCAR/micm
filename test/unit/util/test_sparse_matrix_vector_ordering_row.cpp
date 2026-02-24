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

