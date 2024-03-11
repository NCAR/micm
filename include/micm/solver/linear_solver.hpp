// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <functional>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief Reorders a set of state variables using Diagonal Markowitz algorithm
  /// @param matrix Original matrix non-zero elements
  /// @result Reordered mapping vector (reordered[i] = original[map[i]])
  template<template<class> class MatrixPolicy>
  std::vector<std::size_t> DiagonalMarkowitzReorder(const MatrixPolicy<int>& matrix);

  /// @brief A general-use block-diagonal sparse-matrix linear solver
  ///
  /// The sparsity pattern of each block in the block diagonal matrix is the same.
  template<typename T, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy = LuDecomposition>
  class LinearSolver
  {
   protected:
    // Parameters needed to calculate L (U x) = b
    //
    // The calculation is split into calculation of L y = b where y = U x:
    //
    // y_1 = b_1 / L_11
    // y_i = 1 / L_ii * [ b_i - sum( j = 1...i-1 ){ L_ij * y_j } ]  i = 2...N
    //
    // ... and then U x = y:
    //
    // x_N = y_N / U_NN
    // x_i = 1 / U_ii * [ y_i - sum( j = i+1...N ){ U_ij * x_j } ] i = N-1...1

    // Number of non-zero elements (excluding the diagonal) and the index of the diagonal
    // element for each row in L
    std::vector<std::pair<std::size_t, std::size_t>> nLij_Lii_;
    // Indices of non-zero combinations of L_ij and y_j
    std::vector<std::pair<std::size_t, std::size_t>> Lij_yj_;
    // Number of non-zero elements (exluding the diagonal) and the index of the diagonal
    // element for each row in U (in reverse order)
    std::vector<std::pair<std::size_t, std::size_t>> nUij_Uii_;
    // Indices of non-zero combinations of U_ij and x_j
    std::vector<std::pair<std::size_t, std::size_t>> Uij_xj_;

    LuDecompositionPolicy lu_decomp_;

   public:
    /// @brief default constructor
    LinearSolver(){};

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    LinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value);

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    /// @param create_lu_decomp Function to create an LU Decomposition object that adheres to LuDecompositionPolicy
    LinearSolver(
        const SparseMatrixPolicy<T>& matrix,
        T initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<T>&)> create_lu_decomp);

    /// @brief Decompose the matrix into upper and lower triangular matrices
    void
    Factor(const SparseMatrixPolicy<T>& matrix, SparseMatrixPolicy<T>& lower_matrix, SparseMatrixPolicy<T>& upper_matrix);

    /// @brief Decompose the matrix into upper and lower triangular matrices
    /// @param matrix Matrix to decompose into lower and upper triangular matrices
    /// @param is_singular Flag that is set to true if matrix is singular; false otherwise
    void Factor(
        const SparseMatrixPolicy<T>& matrix,
        SparseMatrixPolicy<T>& lower_matrix,
        SparseMatrixPolicy<T>& upper_matrix,
        bool& is_singular);

    /// @brief Solve for x in Ax = b
    template<template<class> class MatrixPolicy>
    requires(!VectorizableDense<MatrixPolicy<T>> || !VectorizableSparse<SparseMatrixPolicy<T>>) void Solve(
        const MatrixPolicy<T>& b,
        MatrixPolicy<T>& x,
        SparseMatrixPolicy<T>& lower_matrix,
        SparseMatrixPolicy<T>& upper_matrix);
    template<template<class> class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy<T>>&& VectorizableSparse<SparseMatrixPolicy<T>>) void Solve(
        const MatrixPolicy<T>& b,
        MatrixPolicy<T>& x,
        SparseMatrixPolicy<T>& lower_matrix,
        SparseMatrixPolicy<T>& upper_matrix);
  };

}  // namespace micm

#include "linear_solver.inl"