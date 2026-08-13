// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <functional>

namespace micm
{
  /// @brief A general-use block-diagonal sparse-matrix linear solver
  ///
  /// The sparsity pattern of each block in the block diagonal matrix is the same.
  /// The L and U matrices are decomposed in-place over the original A matrix.
  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy = LuDecompositionInPlace<SparseMatrixPolicy>>
  class LinearSolverInPlace
  {
    using SparseMatrix = SparseMatrixPolicy;
    using DenseMatrix = MatrixPolicy;
    template<class U>
    using Vector = typename SparseMatrix::template VectorType<U>;
    template<class U>
    using VectorView = typename SparseMatrix::VectorType<U>::ConstViewType;
   public:
    using DenseMatrixType = MatrixPolicy;
    using SparseMatrixType = SparseMatrixPolicy;
    using LuDecompositionType = LuDecompositionPolicy;
    struct IndexPair
    {
      Index first_;
      Index second_;
    };
    struct Views
    {
      VectorView<Index> nLij_;
      VectorView<IndexPair> Lij_yj_;
      VectorView<IndexPair> nUij_Uii_;
      VectorView<IndexPair> Uij_xj_;

      Views() = default;

      Views(
        const Vector<Index>& nLij,
        const Vector<IndexPair>& Lij_yj,
        const Vector<IndexPair>& nUij_Uii,
        const Vector<IndexPair>& Uij_xj
      ) : nLij_(nLij.GetView()),
          Lij_yj_(Lij_yj.GetView()),
          nUij_Uii_(nUij_Uii.GetView()),
          Uij_xj_(Uij_xj.GetView())
      {
      }
    };
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
    Vector<Index> nLij_;
    // Indices of non-zero combinations of L_ij and y_j
    Vector<IndexPair> Lij_yj_;
    // Number of non-zero elements (exluding the diagonal) and the index of the diagonal
    // element for each row in U (in reverse order)
    Vector<IndexPair> nUij_Uii_;
    // Indices of non-zero combinations of U_ij and x_j
    Vector<IndexPair> Uij_xj_;
    // MICM_LAMBDA compatible views of the index vectors
    Views views_;

    LuDecompositionPolicy lu_decomp_;

   public:
    /// @brief default constructor
    LinearSolverInPlace() = default;

    LinearSolverInPlace(const LinearSolverInPlace&) = delete;
    LinearSolverInPlace& operator=(const LinearSolverInPlace&) = delete;
    LinearSolverInPlace(LinearSolverInPlace&&) noexcept;
    LinearSolverInPlace& operator=(LinearSolverInPlace&&) noexcept;

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
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)>& create_lu_decomp);

    virtual ~LinearSolverInPlace() = default;

    /// @brief Decompose the matrix into upper and lower triangular matrices (matrix will be overwritten)
    /// @param matrix Matrix to decompose in-place into lower and upper triangular matrices
    void Factor(SparseMatrixPolicy& matrix) const;

    /// @brief Solve for x in Ax = b. x should be a copy of b and after Solve finishes x will contain the result
    /// @param x The solution vector
    /// @param LU The LU decomposition of the matrix as a square sparse matrix
    void Solve(MatrixPolicy& x, const SparseMatrixPolicy& lu_matrix) const;
  };

}  // namespace micm

#include "linear_solver_in_place.inl"