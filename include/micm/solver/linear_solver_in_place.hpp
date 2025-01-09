// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cmath>
#include <functional>

namespace micm
{
  /// @brief A general-use block-diagonal sparse-matrix linear solver
  ///
  /// The sparsity pattern of each block in the block diagonal matrix is the same.
  /// The L and U matrices are decomposed in-place over the original A matrix.
  template<class SparseMatrixPolicy, class LuDecompositionPolicy = LuDecompositionInPlace>
  class LinearSolverInPlace
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

    // Number of non-zero elements (excluding the diagonal) for each row in L
    std::vector<std::size_t> nLij_;
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
    LinearSolverInPlace(){};

    LinearSolverInPlace(const LinearSolverInPlace&) = delete;
    LinearSolverInPlace& operator=(const LinearSolverInPlace&) = delete;
    LinearSolverInPlace(LinearSolverInPlace&&) = default;
    LinearSolverInPlace& operator=(LinearSolverInPlace&&) = default;

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    LinearSolverInPlace(const SparseMatrixPolicy& matrix, typename SparseMatrixPolicy::value_type initial_value);

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    /// @param create_lu_decomp Function to create an LU Decomposition object that adheres to LuDecompositionPolicy
    LinearSolverInPlace(
        const SparseMatrixPolicy& matrix,
        typename SparseMatrixPolicy::value_type initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)> create_lu_decomp);

    virtual ~LinearSolverInPlace() = default;

    /// @brief Decompose the matrix into upper and lower triangular matrices (matrix will be overwritten)
    /// @param matrix Matrix to decompose in-place into lower and upper triangular matrices
    void Factor(SparseMatrixPolicy& matrix) const;

    /// @brief Solve for x in Ax = b. x should be a copy of b and after Solve finishes x will contain the result
    /// @param x The solution vector
    /// @param LU The LU decomposition of the matrix as a square sparse matrix
    template<class MatrixPolicy>
      requires(!VectorizableDense<MatrixPolicy> || !VectorizableSparse<SparseMatrixPolicy>)
    void Solve(MatrixPolicy& x, const SparseMatrixPolicy& lu_matrix) const;
    template<class MatrixPolicy>
      requires(VectorizableDense<MatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    void Solve(MatrixPolicy& x, const SparseMatrixPolicy& lu_matrix) const;
  };

}  // namespace micm

#include "linear_solver_in_place.inl"