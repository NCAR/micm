// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "solver_builder.hpp"

#include <micm/process/cuda_process_set.hpp>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/util/cuda_dense_matrix.hpp>
#include <micm/util/cuda_sparse_matrix.hpp>

namespace micm
{
  /// @brief Builder of CUDA-based general solvers
  /// @tparam SolverParametersPolicy Policy for the ODE solver
  /// @tparam T Primary data type
  /// @tparam L Vector size
  ///
  /// GPU solvers only work with vector-ordered matrices
  template<class SolverParametersPolicy, class T = double, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  using CudaSolverBuilder = SolverBuilder<SolverParametersPolicy, CudaVectorMatrix<T, L>, CudaSparseMatrix<T, SparseMatrixVectorOrdering<L>>, CudaProcessSet, CudaLinearSolver<CudaSparseMatrix<T, SparseMatrixVectorOrdering<L>>, CudaLuDecomposition>>;
} // namespace micm