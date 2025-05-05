// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix following the Doolittle algorithm
  ///
  /// The LU decomposition uses the Doolittle algorithm following the
  /// naming used here: https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
  ///
  /// The sudo-code for the corresponding dense matrix algorithm for matrix A
  /// and lower (upper) triangular matrix L(U) would be:
  ///
  /// for i = 0...n-1                 // Outer loop over rows (columns) for upper (lower) triangular matrix
  ///   for k = i...n-1               // Middle loop over columns for upper triangular matrix
  ///     sum = 0
  ///     for j = 0...i-1             // Inner loop over columns (rows) for lower (upper) triangular matrix
  ///       sum += L[i][j] * U[j][k]
  ///     U[i][k] = A[i][k] - sum
  ///   L[i][i] = 1                   // Lower triangular matrix is 1 along the diagonal
  ///   for k = i+1...n-1             // Middle loop over rows for lower triangular matrix
  ///     sum = 0
  ///     for j = 0...i-1             // Inner loop over columns (rows) for lower (upper) triangular matrix
  ///       sum += L[k][j] * U[j][i];
  ///     L[k][i] = (A[k][i] - sum) / U[i][i]
  ///
  /// For the sparse matrix algorithm, the indices of non-zero terms are stored in
  /// several arrays during construction. These arrays are iterated through during
  /// calls to Decompose to do the actual decomposition.
  /// Our LU Decomposition only assigns the values of the jacobian to the LU matrices
  /// when the *jacobian* is nonzero. However, the sparsity pattern of the jacobian doesn't
  /// necessarily match that of the LU matrices. There can be more nonzero elements in the LU matrices
  /// than in the jacobian. When this happens, we still need to assign the value of the jacobian matrix
  /// to the LU matrix. This value is implicitly zero when the sparsity pattern differs. The Fill values
  /// here do this implicit assignment
  /// More detail in this issue: https://github.com/NCAR/micm/issues/625
  class LuDecompositionDoolittle
  {
   protected:
    /// number of elements in the middle (k) loops for lower and upper triangular matrices, respectively,
    /// for each iteration of the outer (i) loop
    std::vector<std::pair<std::size_t, std::size_t>> niLU_;
    /// True when A[i][k] is non-zero for each iteration of the middle (k) loop for the upper
    /// triangular matrix; False otherwise. Used data type char instead of bool because vector<bool> representation
    /// does not support easy retrieval of memory address using data() function.
    std::vector<char> do_aik_;
    /// Index in A.data_ for A[i][k] for each iteration of the middle (k) loop for the upper
    /// triangular matrix when A[i][k] is non-zero
    std::vector<std::size_t> aik_;
    /// Index in U.data_ for U[i][k] for each iteration of the middle (k) loop for the upper
    /// triangular matrix when U[i][k] is non-zero, and the corresponding number of elements
    /// in the inner (j) loop
    std::vector<std::pair<std::size_t, std::size_t>> uik_nkj_;
    /// Index in L.data_ for L[i][j], and in U.data_ for U[j][k] in the upper inner (j) loop
    /// when L[i][j] and U[j][k] are both non-zero.
    std::vector<std::pair<std::size_t, std::size_t>> lij_ujk_;
    /// True when A[k][i] is non-zero for each iteration of the middle (k) loop for the lower
    /// triangular matrix; False otherwise. Used data type char instead of bool because vector<bool> representation
    /// does not suppor easy retrieval of memory address using data() function.
    std::vector<char> do_aki_;
    /// Index in A.data_ for A[k][i] for each iteration of the middle (k) loop for the lower
    /// triangular matrix when A[k][i] is non-zero.
    std::vector<std::size_t> aki_;
    /// Index in L.data_ for L[k][i] for each iteration of the middle (k) loop for the lower
    /// triangular matrix when L[k][i] is non-zero, and the corresponding number of elements
    /// in the inner (j) loop
    std::vector<std::pair<std::size_t, std::size_t>> lki_nkj_;
    /// Index in L.data_ for L[k][j], and in U.data_ for U[j][i] in the lower inner (j) loop
    /// when L[k][j] and U[j][i] are both non-zero.
    std::vector<std::pair<std::size_t, std::size_t>> lkj_uji_;
    /// Index in U.data_ for U[i][i] for each interation in the middle (k) loop for the lower
    /// triangular matrix when L[k][i] is non-zero
    std::vector<std::size_t> uii_;

   public:
    /// @brief default constructor
    LuDecompositionDoolittle();

    LuDecompositionDoolittle(const LuDecompositionDoolittle&) = delete;
    LuDecompositionDoolittle& operator=(const LuDecompositionDoolittle&) = delete;

    LuDecompositionDoolittle(LuDecompositionDoolittle&& other) = default;
    LuDecompositionDoolittle& operator=(LuDecompositionDoolittle&&) = default;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    LuDecompositionDoolittle(const SparseMatrixPolicy& matrix);

    ~LuDecompositionDoolittle() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static LuDecompositionDoolittle Create(const SparseMatrixPolicy& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return L and U Sparse matrices
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static std::pair<LMatrixPolicy, UMatrixPolicy> GetLUMatrices(
        const SparseMatrixPolicy& A,
        typename SparseMatrixPolicy::value_type initial_value);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param L The lower triangular matrix created by decomposition
    /// @param U The upper triangular matrix created by decomposition
    template<class SparseMatrixPolicy>
      requires(!VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const;
    template<class SparseMatrixPolicy>
      requires(VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const;

   protected:
    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_doolittle.inl"
