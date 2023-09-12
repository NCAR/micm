// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<std::size_t L, template<class> class SparseMatrixPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy>::JitLinearSolver(
      std::shared_ptr<JitCompiler> compiler,
      const SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
      double initial_value)
      : LinearSolver<double, SparseMatrixPolicy>(matrix, initial_value),
        compiler_(compiler)
  {
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy>
  inline JitLinearSolver<L, SparseMatrixPolicy>::~JitLinearSolver()
  {
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy>::Factor(SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix)
  {
    LinearSolver<double, SparseMatrixPolicy>::Factor(matrix);
  }

  template<std::size_t L, template<class> class SparseMatrixPolicy>
  template<template<class> class MatrixPolicy>
  inline void JitLinearSolver<L, SparseMatrixPolicy>::Solve(const MatrixPolicy<double>& b, MatrixPolicy<double>& x)
  {
    LinearSolver<double, SparseMatrixPolicy>::template Solve<MatrixPolicy>(b, x);
  }

}  // namespace micm