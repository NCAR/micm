// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../util/test_sparse_matrix_policy.hpp"

#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

using KokkosVectorOrdering1 = micm::SparseMatrixVectorOrderingCompressedSparseRow<1>;
using KokkosVectorOrdering2 = micm::SparseMatrixVectorOrderingCompressedSparseRow<2>;
using KokkosVectorOrdering3 = micm::SparseMatrixVectorOrderingCompressedSparseRow<3>;
using KokkosVectorOrdering4 = micm::SparseMatrixVectorOrderingCompressedSparseRow<4>;

TEST(KokkosSparseMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testZeroMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering2>();
}

TEST(KokkosSparseMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testConstZeroMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering3>();
}

TEST(KokkosSparseMatrix, SetScalar)
{
  testSetScalar<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testSetScalar<micm::KokkosSparseMatrix, KokkosVectorOrdering3>();
}

TEST(KokkosSparseMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosVectorOrdering3>();
  testAddToDiagonal<micm::KokkosSparseMatrix, KokkosVectorOrdering4>();
}

TEST(KokkosSparseMatrix, SingleBlockMatrix)
{
  testSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering4>();
}

TEST(KokkosSparseMatrix, ConstSingleBlockMatrix)
{
  testConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testConstSingleBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering2>();
}

TEST(KokkosSparseMatrix, MultiBlockMatrix)
{
  testMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testMultiBlockMatrix<micm::KokkosSparseMatrix, KokkosVectorOrdering2>();
}

TEST(KokkosSparseMatrix, Print)
{
  testPrint<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testPrint<micm::KokkosSparseMatrix, KokkosVectorOrdering2>();
  testPrint<micm::KokkosSparseMatrix, KokkosVectorOrdering3>();
  testPrint<micm::KokkosSparseMatrix, KokkosVectorOrdering4>();
}

TEST(KokkosSparseMatrix, PrintNonZero)
{
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosVectorOrdering1>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosVectorOrdering2>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosVectorOrdering3>();
  testPrintNonZero<micm::KokkosSparseMatrix, KokkosVectorOrdering4>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
