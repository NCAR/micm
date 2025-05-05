// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix following the algorithm from the MOZART model
  ///
  /// This LU decomposition uses the algorithm from the MOZART chemistry preprocessor at:
  /// https://github.com/ESCOMP/CHEM_PREPROCESSOR/blob/f923081508f4264e61fcef2753a9898e55d1598e/src/cam_chempp/make_lu_fac.f#L127-L214
  ///
  /// The MOZART function overwrote the A matrix with the L and U matrices. The pseudo-code
  /// in C++ for the corresponding dense matrix algorithm for matrix A (inline change) would be:
  ///
  /// for i = 0...n-1                     // Outer loop over columns of the sparse matrix A
  ///   for j = i+1...n-1                 // Multiply column below diagonal
  ///     A[j][i] = A[j][i] / A[i][i]
  ///   for k = i+1...n-1                 // Modify sub-matrix
  ///     for j = i+1...n-1
  ///       A[j][k] = A[j][k] – A[j][i] * A[i][k]
  ///
  /// The pseudo-code in C++ for the corresponding dense matrix algorithm for matrix A
  /// and separate lower (upper) triangular matrix L(U) would be:
  ///
  /// for i = 0...n-1                     // Initialize U and L matrices to the A values
  ///   for j = 0...i                     // Initialize U matrix including diagonal
  ///     U[j][i] = A[j][i]
  ///   L[i][i] = 1                       // Lower triangular matrix is 1 along the diagonal
  ///   for j = i+1...n-1                 // Initialize L matrix excluding diagonal
  ///     L[j][i] = A[j][i]
  /// for i = 0...n-1
  ///   for j = i+1...n-1                 // Multiply column below diagonal
  ///     L[j][i] = L[j][i] / U[i][i]
  ///   for k = i+1...n-1                 // Modify sub-matrix
  ///     for j = i+1...k
  ///       U[j][k] = U[j][k] - L[j][i] * U[i][k]
  ///     for j = k+1...n-1
  ///       L[j][k] = L[j][k] - L[j][i] * U[i][k]
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
  class LuDecompositionMozart
  {
   protected:
    /// Index in L.data_ for all diagonal elements, and number of iterations of the middle (j) loops
    /// used to set the initial value for the L and U matrices
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> lii_nuji_nlji_;
    /// Index in U.data_ and A.data_ for U[j][i] and A[j][i] for each iteration of the inner (j) loop
    /// used to set the initial value for the U matrix
    std::vector<std::pair<std::size_t, std::size_t>> uji_aji_;
    /// Index in L.data_ and A.data_ for L[j][i] and A[j][i] for each iteration of the inner (j) loop
    /// used to set the initial value for the L matrix
    std::vector<std::pair<std::size_t, std::size_t>> lji_aji_;
    /// Index in U.data_ for each non-zero element in U that is zero in A
    std::vector<std::size_t> fill_uji_;
    /// Index in L.data_ for each non-zero element in L that is zero in A
    std::vector<std::size_t> fill_lji_;
    /// Index in U.data_ for U[i][i] and the number of elements in the middle (j) and (k) loops
    /// for each iteration of the outer (i) loop
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> uii_nj_nk_;
    /// Index in L.data_ for L[j][i] for each iteration of the middle (j) loop
    /// for the lower triangular matrix
    std::vector<std::size_t> lji_;
    /// Number of elements in the inner (j) loops for each iteration of the middle (k) loop for the
    /// upper and lower triangular matrices, and the index in U.data_ for U[j][k] for each iteration
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> nujk_nljk_uik_;
    /// Index in U.data_ for U[j][k] and in L.data_ for L[j][i]
    /// for each iteration of the inner (j) loop for the upper triangular matrix
    std::vector<std::pair<std::size_t, std::size_t>> ujk_lji_;
    /// Index in L.data_ for L[j][k] and in L.data_ for L[j][i]
    /// for each iteration of the inner (j) loop for the lower triangular matrix
    std::vector<std::pair<std::size_t, std::size_t>> ljk_lji_;

   public:
    /// @brief default constructor
    LuDecompositionMozart();

    LuDecompositionMozart(const LuDecompositionMozart&) = delete;
    LuDecompositionMozart& operator=(const LuDecompositionMozart&) = delete;

    LuDecompositionMozart(LuDecompositionMozart&& other) = default;
    LuDecompositionMozart& operator=(LuDecompositionMozart&&) = default;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    LuDecompositionMozart(const SparseMatrixPolicy& matrix);

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    LuDecompositionMozart(const SparseMatrixPolicy& matrix);

    ~LuDecompositionMozart() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static LuDecompositionMozart Create(const SparseMatrixPolicy& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return L and U Sparse matrices
    template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
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

   private:
    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_mozart.inl"
