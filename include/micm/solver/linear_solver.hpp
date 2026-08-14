// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/solver/linear_solver_in_place.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <functional>

namespace micm
{

  /// @brief Concept for in-place linear solver algorithms
  template<class T, class DenseMatrixPolicy, class SparseMatrixPolicy>
  concept LinearSolverInPlaceConcept = requires(T t) {
    { t.Factor(std::declval<SparseMatrixPolicy&>()) } -> std::same_as<void>;
    { t.Solve(std::declval<DenseMatrixPolicy&>(), std::declval<SparseMatrixPolicy>()) } -> std::same_as<void>;
  };
  static_assert(
      LinearSolverInPlaceConcept<
          LinearSolverInPlace<StandardDenseMatrix, StandardSparseMatrix>,
          StandardDenseMatrix,
          StandardSparseMatrix>,
      "LinearSolverInPlace does not meet the LinearSolverInPlaceConcept requirements");
  static_assert(
      LinearSolverInPlaceConcept<
          LinearSolverInPlace<
              VectorMatrix<Real, 1>,
              SparseMatrix<Real, SparseMatrixVectorOrderingCompressedSparseRow<1>>,
              LuDecompositionMozartInPlace<SparseMatrix<Real, SparseMatrixVectorOrderingCompressedSparseRow<1>>>>,
          VectorMatrix<Real, 1>,
          SparseMatrix<Real, SparseMatrixVectorOrderingCompressedSparseRow<1>>>,
      "LinearSolverInPlace for vector matrices does not meet the LinearSolverInPlaceConcept requirements");

  /// @brief Reorders a set of state variables using Diagonal Markowitz algorithm
  /// @param matrix Original matrix non-zero elements
  /// @result Reordered mapping vector (reordered[i] = original[map[i]])
  template<template<class> class MatrixPolicy>
  std::vector<Index> DiagonalMarkowitzReorder(const MatrixPolicy<int>& matrix);

  /// @brief A general-use block-diagonal sparse-matrix linear solver
  ///
  /// The sparsity pattern of each block in the block diagonal matrix is the same.
  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy = LuDecomposition<SparseMatrixPolicy>>
  class LinearSolver
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
    struct Views
    {
      VectorView<IndexPair> nLij_Lii_;
      VectorView<IndexPair> Lij_yj_;
      VectorView<IndexPair> nUij_Uii_;
      VectorView<IndexPair> Uij_xj_;

      Views() = default;

      Views(
        const Vector<IndexPair>& nLij_Lii,
        const Vector<IndexPair>& Lij_yj,
        const Vector<IndexPair>& nUij_Uii,
        const Vector<IndexPair>& Uij_xj
      ) : nLij_Lii_(nLij_Lii.GetView()),
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

    // Number of non-zero elements (excluding the diagonal) and the index of the diagonal
    // element for each row in L
    Vector<IndexPair> nLij_Lii_;
    // Indices of non-zero combinations of L_ij and y_j
    Vector<IndexPair> Lij_yj_;
    // Number of non-zero elements (exluding the diagonal) and the index of the diagonal
    // element for each row in U (in reverse order)
    Vector<IndexPair> nUij_Uii_;
    // Indices of non-zero combinations of U_ij and x_j
    Vector<IndexPair> Uij_xj_;
    // MICM_LAMBDA compatible views of index vectors
    Views views_;

    LuDecompositionPolicy lu_decomp_;

   public:
    /// @brief default constructor
    LinearSolver() = default;

    LinearSolver(const LinearSolver&) = delete;
    LinearSolver& operator=(const LinearSolver&) = delete;
    LinearSolver(LinearSolver&&) noexcept;
    LinearSolver& operator=(LinearSolver&&) noexcept;

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    LinearSolver(const SparseMatrixPolicy& matrix, typename SparseMatrixPolicy::value_type initial_value);

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    /// @param initial_value Initial value for matrix elements
    /// @param create_lu_decomp Function to create an LU Decomposition object that adheres to LuDecompositionPolicy
    LinearSolver(
        const SparseMatrixPolicy& matrix,
        typename SparseMatrixPolicy::value_type initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)>& create_lu_decomp);

    virtual ~LinearSolver() = default;

    /// @brief Decompose the matrix into upper and lower triangular matrices
    /// @param matrix Matrix to decompose into lower and upper triangular matrices
    void Factor(const SparseMatrixPolicy& matrix, SparseMatrixPolicy& lower_matrix, SparseMatrixPolicy& upper_matrix) const;

    /// @brief Solve for x in Ax = b. x should be a copy of b and after Solve finishes x will contain the result
    void Solve(MatrixPolicy& x, const SparseMatrixPolicy& lower_matrix, const SparseMatrixPolicy& upper_matrix) const;
  };

}  // namespace micm

#include "linear_solver.inl"