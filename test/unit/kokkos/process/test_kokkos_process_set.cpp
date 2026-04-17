// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../../process/test_process_set_policy.hpp"

#include <micm/kokkos/process/kokkos_process_set.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>

#include <gtest/gtest.h>

using KokkosDenseMatrix1 = micm::KokkosDenseMatrix<double, 1>;
using KokkosDenseMatrix2 = micm::KokkosDenseMatrix<double, 2>;
using KokkosDenseMatrix3 = micm::KokkosDenseMatrix<double, 3>;
using KokkosDenseMatrix4 = micm::KokkosDenseMatrix<double, 4>;

using KokkosSparseMatrix1 = micm::KokkosSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using KokkosSparseMatrix2 = micm::KokkosSparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using KokkosSparseMatrix3 = micm::KokkosSparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using KokkosSparseMatrix4 = micm::KokkosSparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(KokkosProcessSet, VectorMatrix)
{
  testProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1, micm::KokkosProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1>>();
  testProcessSet<KokkosDenseMatrix2, KokkosSparseMatrix2, micm::KokkosProcessSet<KokkosDenseMatrix2, KokkosSparseMatrix2>>();
  testProcessSet<KokkosDenseMatrix3, KokkosSparseMatrix3, micm::KokkosProcessSet<KokkosDenseMatrix3, KokkosSparseMatrix3>>();
  testProcessSet<KokkosDenseMatrix4, KokkosSparseMatrix4, micm::KokkosProcessSet<KokkosDenseMatrix4, KokkosSparseMatrix4>>();
}

TEST(KokkosRandomProcessSet, VectorMatrix)
{
  testRandomSystem<KokkosDenseMatrix1, KokkosSparseMatrix1, micm::KokkosProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1>>(
      200, 50, 40);
  testRandomSystem<KokkosDenseMatrix1, KokkosSparseMatrix1, micm::KokkosProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1>>(
      300, 30, 20);
  testRandomSystem<KokkosDenseMatrix1, KokkosSparseMatrix1, micm::KokkosProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1>>(
      400, 100, 80);
}

TEST(KokkosProcessSetAlgebraicVariables, VectorMatrix)
{
  testAlgebraicMasking<
      KokkosDenseMatrix1,
      KokkosSparseMatrix1,
      micm::KokkosProcessSet<KokkosDenseMatrix1, KokkosSparseMatrix1>>();
  testAlgebraicMasking<
      KokkosDenseMatrix2,
      KokkosSparseMatrix2,
      micm::KokkosProcessSet<KokkosDenseMatrix2, KokkosSparseMatrix2>>();
  testAlgebraicMasking<
      KokkosDenseMatrix3,
      KokkosSparseMatrix3,
      micm::KokkosProcessSet<KokkosDenseMatrix3, KokkosSparseMatrix3>>();
  testAlgebraicMasking<
      KokkosDenseMatrix4,
      KokkosSparseMatrix4,
      micm::KokkosProcessSet<KokkosDenseMatrix4, KokkosSparseMatrix4>>();
}
