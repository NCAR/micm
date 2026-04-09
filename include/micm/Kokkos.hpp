// Copyright (C) 2024-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/kokkos/process/kokkos_process_set.hpp>
#include <micm/kokkos/util/kokkos_util.hpp>

namespace micm
{
  using KokkosDenseMatrixVector = KokkosDenseMatrix<double, MICM_DEFAULT_VECTOR_SIZE>;
  using KokkosSparseMatrixVector = KokkosSparseMatrix<double, SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>>;

  using KokkosProcessSetVector = KokkosProcessSet<KokkosDenseMatrixVector, KokkosSparseMatrixVector>;

}  // namespace micm
