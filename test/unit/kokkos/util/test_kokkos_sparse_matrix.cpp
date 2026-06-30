// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../util/test_sparse_matrix_policy.hpp"

#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

using KokkosOrdering1 = micm::SparseMatrixVectorOrderingCompressedSparseRow<1>;
using KokkosOrdering2 = micm::SparseMatrixVectorOrderingCompressedSparseRow<2>;
using KokkosOrdering3 = micm::SparseMatrixVectorOrderingCompressedSparseRow<3>;
using KokkosOrdering4 = micm::SparseMatrixVectorOrderingCompressedSparseRow<4>;

// Core matrix tests
TEST(KokkosSparseMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testConstZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering3>();
}

TEST(KokkosSparseMatrix, SetScalar)
{
  testSetScalar<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testSetScalar<micm::KokkosSparseMatrix, KokkosOrdering3>();
}

TEST(KokkosSparseMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, SingleBlockMatrix)
{
  testSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, ConstSingleBlockMatrix)
{
  testConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, MultiBlockMatrix)
{
  testMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, Print)
{
  testPrint<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testPrint<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testPrint<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testPrint<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, PrintNonZero)
{
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// BlockFunction infrastructure tests
TEST(KokkosSparseMatrix, ArrayFunction)
{
  testArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultiMatrixArrayFunction)
{
  testMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultipleTemporaries)
{
  testMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, BlockViewReuse)
{
  testBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, FunctionReusability)
{
  testFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, TwoSparseMatricesDifferentStructure)
{
  testTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Sparse + Dense/Vector matrix interaction tests
TEST(KokkosSparseMatrix, SparseAndDenseMatrixFunction)
{
  testSparseAndDenseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1, micm::Matrix>();
}

TEST(KokkosSparseMatrix, SparseAndVectorMatrixFunction)
{
  testSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1, 1>();
  testSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2, 2>();
  testSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3, 3>();
  testSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4, 4>();
}

// Ordering compatibility tests
TEST(KokkosSparseMatrix, IncompatibleOrdering)
{
  testIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, micm::Matrix>();
  testIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering3, micm::Matrix>();
  testIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering4, micm::Matrix>();
}

TEST(KokkosSparseMatrix, IncompatibleVectorOrdering)
{
  testIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, 1>();
  testIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, 3>();
  testIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering4, 2>();
}

TEST(KokkosSparseMatrix, IncompatibleSparseOrdering)
{
  testIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering1, KokkosOrdering2>();
  testIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, KokkosOrdering3>();
  testIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering3, KokkosOrdering4>();
}

// Block dimension mismatch tests
TEST(KokkosSparseMatrix, MismatchedBlockDimensions)
{
  testMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MismatchedElementDimensions)
{
  testMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, WrongMatrixDimensions)
{
  testWrongMatrixDimensions<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testWrongMatrixDimensions<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testWrongMatrixDimensions<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testWrongMatrixDimensions<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, ConstSparseMatrixFunction)
{
  testConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, EmptySparseMatrixFunction)
{
  testEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Different blocks from creation tests
TEST(KokkosSparseMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Mismatched blocks at invocation tests
TEST(KokkosSparseMatrix, MismatchedBlocksAtInvocation)
{
  testMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, WrongStructureAtInvocation)
{
  testWrongStructureAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testWrongStructureAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testWrongStructureAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testWrongStructureAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Vector support tests
TEST(KokkosSparseMatrix, VectorInSparseMatrixFunction)
{
  testVectorInSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorTooSmall)
{
  testVectorTooSmall<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorTooLarge)
{
  testVectorTooLarge<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, EmptyVectorNonEmptySparseMatrix)
{
  testEmptyVectorNonEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, NonEmptyVectorEmptySparseMatrix)
{
  testNonEmptyVectorEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, EmptyVectorEmptySparseMatrix)
{
  testEmptyVectorEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleVectorsDifferentSizes)
{
  testMultipleVectorsDifferentSizes<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleVectorsSameSize)
{
  testMultipleVectorsSameSize<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesOneVector)
{
  testMultipleSparseMatricesOneVector<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  testMultipleSparseMatricesDifferentBlocksVector<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  testVectorSizeMatchesOneSparseMatrixOnly<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, ConstVectorSparse)
{
  testConstVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MutableVectorSparse)
{
  testMutableVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, FunctionReusabilityWithVectorsSparse)
{
  testFunctionReusabilityWithVectorsSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  testFunctionInvocationWithWrongSizedVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, ArraySupportSparse)
{
  testArraySupportSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MixedVectorBlockViewBlockVariable)
{
  testMixedVectorBlockViewBlockVariable<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, IntegerVectorSparse)
{
  testIntegerVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, FunctionWithConstSignatureSparse)
{
  testFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, GetBlockViewByVectorIndex)
{
  testGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering1>();
  testGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering2>();
  testGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering3>();
  testGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
