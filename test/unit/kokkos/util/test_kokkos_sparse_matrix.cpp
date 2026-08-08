// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../util/test_sparse_matrix_policy.hpp"

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

using KokkosOrdering1 = micm::SparseMatrixVectorOrderingCompressedSparseRow<1>;
using KokkosOrdering2 = micm::SparseMatrixVectorOrderingCompressedSparseRow<2>;
using KokkosOrdering3 = micm::SparseMatrixVectorOrderingCompressedSparseRow<3>;
using KokkosOrdering4 = micm::SparseMatrixVectorOrderingCompressedSparseRow<4>;

template<class T>
using KokkosDense1 = micm::KokkosDenseMatrix<T,1>;
template<class T>
using KokkosDense2 = micm::KokkosDenseMatrix<T,2>;
template<class T>
using KokkosDense3 = micm::KokkosDenseMatrix<T,3>;
template<class T>
using KokkosDense4 = micm::KokkosDenseMatrix<T,4>;


// Core matrix tests
TEST(KokkosSparseMatrix, ZeroMatrix)
{
  TestZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, ConstZeroMatrix)
{
  TestConstZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestConstZeroMatrix<micm::KokkosSparseMatrix, KokkosOrdering3>();
}

TEST(KokkosSparseMatrix, SetScalar)
{
  TestSetScalar<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestSetScalar<micm::KokkosSparseMatrix, KokkosOrdering3>();
}

TEST(KokkosSparseMatrix, AddToDiagonal)
{
  TestAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestAddToDiagonal<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, SingleBlockMatrix)
{
  TestSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, ConstSingleBlockMatrix)
{
  TestConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, MultiBlockMatrix)
{
  TestMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosOrdering2>();
}

TEST(KokkosSparseMatrix, Print)
{
  TestPrint<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestPrint<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestPrint<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestPrint<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, PrintNonZero)
{
  TestPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestPrintNonZero<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// BlockFunction infrastructure tests
TEST(KokkosSparseMatrix, ArrayFunction)
{
  TestArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultiMatrixArrayFunction)
{
  TestMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMultiMatrixArrayFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultipleTemporaries)
{
  TestMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMultipleTemporaries<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, BlockViewReuse)
{
  TestBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestBlockViewReuse<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, FunctionReusability)
{
  TestFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestFunctionReusability<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, TwoSparseMatricesDifferentStructure)
{
  TestTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestTwoSparseMatricesDifferentStructure<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, SparseAndDenseMatrixFunction)
{
  TestSparseAndDenseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1, KokkosDense1>();
  TestSparseAndDenseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2, KokkosDense2>();
  TestSparseAndDenseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3, KokkosDense3>();
  TestSparseAndDenseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4, KokkosDense4>();
}

TEST(KokkosSparseMatrix, SparseAndVectorMatrixFunction)
{
  TestSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1, KokkosDense1, 1>();
  TestSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2, KokkosDense2, 2>();
  TestSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3, KokkosDense3, 3>();
  TestSparseAndVectorMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4, KokkosDense4, 4>();
}

#if 0

// These tests check for exceptions, but incompatible ordering for Kokkos matrix functions
// is a compile-time check.

// Ordering compatibility tests
TEST(KokkosSparseMatrix, IncompatibleOrdering)
{
  TestIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, KokkosDense1>();
  TestIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering3, KokkosDense4>();
  TestIncompatibleOrdering<micm::KokkosSparseMatrix, KokkosOrdering4, KokkosDense2>();
}

TEST(KokkosSparseMatrix, IncompatibleVectorOrdering)
{
  TestIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, 1>();
  TestIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, 3>();
  TestIncompatibleVectorOrdering<micm::KokkosSparseMatrix, KokkosOrdering4, 2>();
}

TEST(KokkosSparseMatrix, IncompatibleSparseOrdering)
{
  TestIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering1, KokkosOrdering2>();
  TestIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering2, KokkosOrdering3>();
  TestIncompatibleSparseOrdering<micm::KokkosSparseMatrix, KokkosOrdering3, KokkosOrdering4>();
}

#endif

// Block dimension mismatch tests
TEST(KokkosSparseMatrix, MismatchedBlockDimensions)
{
  TestMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMismatchedBlockDimensions<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MismatchedElementDimensions)
{
  TestMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMismatchedElementDimensions<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, ConstSparseMatrixFunction)
{
  TestConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestConstSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, EmptySparseMatrixFunction)
{
  TestEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestEmptySparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Different blocks from creation tests
TEST(KokkosSparseMatrix, MultipleSparseMatricesDifferentBlocksFromCreation)
{
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMultipleSparseMatricesDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, SparseMatrixVectorDifferentBlocksFromCreation)
{
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestSparseMatrixVectorDifferentBlocksFromCreation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Mismatched blocks at invocation tests
TEST(KokkosSparseMatrix, MismatchedBlocksAtInvocation)
{
  TestMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesMismatchedBlocksAtInvocation)
{
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestMultipleSparseMatricesMismatchedBlocksAtInvocation<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

// Vector support tests
TEST(KokkosSparseMatrix, VectorInSparseMatrixFunction)
{
  TestVectorInSparseMatrixFunction<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorTooSmall)
{
  TestVectorTooSmall<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorTooLarge)
{
  TestVectorTooLarge<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, EmptyVectorNonEmptySparseMatrix)
{
  TestEmptyVectorNonEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, NonEmptyVectorEmptySparseMatrix)
{
  TestNonEmptyVectorEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, EmptyVectorEmptySparseMatrix)
{
  TestEmptyVectorEmptySparseMatrix<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleVectorsDifferentSizes)
{
  TestMultipleVectorsDifferentSizes<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleVectorsSameSize)
{
  TestMultipleVectorsSameSize<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesOneVector)
{
  TestMultipleSparseMatricesOneVector<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MultipleSparseMatricesDifferentBlocksVector)
{
  TestMultipleSparseMatricesDifferentBlocksVector<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, VectorSizeMatchesOneSparseMatrixOnly)
{
  TestVectorSizeMatchesOneSparseMatrixOnly<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, ConstVectorSparse)
{
  TestConstVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MutableVectorSparse)
{
  TestMutableVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, FunctionReusabilityWithVectorsSparse)
{
  TestFunctionReusabilityWithVectorsSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, FunctionInvocationWithWrongSizedVectorSparse)
{
  TestFunctionInvocationWithWrongSizedVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, MixedVectorBlockViewBlockVariable)
{
  TestMixedVectorBlockViewBlockVariable<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

TEST(KokkosSparseMatrix, IntegerVectorSparse)
{
  TestIntegerVectorSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
}

#ifndef KOKKOS_ENABLE_CUDA
// This test attemps to wrap the function in std::function which isn't allowed
// with CUDA

TEST(KokkosSparseMatrix, FunctionWithConstSignatureSparse)
{
  TestFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestFunctionWithConstSignatureSparse<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

#endif

TEST(KokkosSparseMatrix, GetBlockViewByVectorIndex)
{
  TestGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering1>();
  TestGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering2>();
  TestGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering3>();
  TestGetBlockViewByVectorIndex<micm::KokkosSparseMatrix, KokkosOrdering4>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
