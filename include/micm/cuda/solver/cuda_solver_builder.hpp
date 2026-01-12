/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/cuda/process/cuda_process_set.hpp>
#include <micm/cuda/solver/cuda_linear_solver_in_place.hpp>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/solver/cuda_state.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{
  /// @brief Builder of CUDA-based general solvers
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam L Vector size
  ///
  /// GPU solvers only work with vector-ordered matrices
  template<class SolverParametersPolicy, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using CudaSolverBuilderInPlace = SolverBuilder<
      SolverParametersPolicy,
      CudaDenseMatrix<double, L>,
      CudaSparseMatrix<double, SparseMatrixVectorOrdering<L>>,
      CudaProcessSet,
      CudaLuDecompositionMozartInPlace,
      CudaLinearSolverInPlace<CudaSparseMatrix<double, SparseMatrixVectorOrdering<L>>, CudaLuDecompositionMozartInPlace>,
      CudaState<
          CudaDenseMatrix<double, L>,
          CudaSparseMatrix<double, SparseMatrixVectorOrdering<L>>,
          CudaLuDecompositionMozartInPlace>>;
}  // namespace micm