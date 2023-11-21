// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <functional>
#include <micm/solver/cuda_linear_solver.cuh>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  template<typename T, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy = CudaLuDecomposition>
  class CudaLinearSolver : public LinearSolver<T, SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    // constructor
    CudaLinearSolver(){};

    CudaLinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value)
        : LinearSolver<T, SparseMatrixPolicy, LuDecompositionPolicy>(matrix, initial_value){};

    CudaLinearSolver(
        const SparseMatrixPolicy<T>& matrix,
        T initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<T>&)> create_lu_decomp)
        : linearSolver(matrix, initial_value, create_lu_decomp);

    template<template<class> class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy<T>> || VectorizableSparse<SparseMatrixPolicy<T>>) std::chrono::nanoseconds Solve(
        const MatrixPolicy<T>& b,
        MatrixPolicy<T>& x,
        SparseMatrixPolicy<T>& lower_matrix,
        SparseMatrixPolicy<T>& upper_matrix)
    {
      CudaLinearSolverParam linearSolver;
      CudaSparseMatrixParam sparseMatrix;
      CudaMatrixParam denseMatrix;
      linearSolver.nLij_Lii_ = this->nLij_Lii_.data();
      linearSolver.nLij_Lii_size_ = this->nLij_Lii_.size();
      linearSolver.Lij_yj_ = this->Lij_yj_.data();
      linearSolver.Lij_yj_size_ = this->Lij_yj_.size();
      linearSolver.nUij_Uii_ = this->nUij_Uii_.data();
      linearSolver.nUij_Uii_size_ = this->nUij_Uii_.size();
      linearSolver.Uij_xj_ = this->Uij_xj_.data();
      linearSolver.Uij_xj_size_ = this->Uij_xj_.size();
      sparseMatrix.lower_matrix_ = lower_matrix.AsVector().data();
      sparseMatrix.lower_matrix_size_ = lower_matrix.AsVector().size();
      sparseMatrix.upper_matrix_ = upper_matrix.AsVector().data();
      sparseMatrix.upper_matrix_size_ = upper_matrix.AsVector().size();
      denseMatrix.b_ = b.AsVector().data();
      denseMatrix.x_ = x.AsVector().data();
      denseMatrix.b_size_ = b.AsVector().size();
      denseMatrix.x_size_ = x.AsVector().size();
      denseMatrix.n_grids_ = b.size();  // number of grids
      denseMatrix.b_column_counts_ = b[0].size();
      denseMatrix.x_column_counts_ = x[0].size();
      // calling kernel driver
      return micm::cuda::SolveKernelDriver(linearSolver, sparseMatrix, denseMatrix);
    };
  };
}  // namespace micm
