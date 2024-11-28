// Copyright (C) 2023-2024 National Center for Atmospheric Research
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
  /// in C++ for the corresponding sparse matrix algorithm for matrix A (inline change) would be: 
  ///
  /// for i = 0...n-1                     // Outer loop over columns of the sparse matrix A
  ///   A[i][i] = 1 / A[i][i]             // Form diagonal inverse
  ///   for j = i+1...n-1                 // Multiply column below diagonal
  ///     A[j][i] = A[j][i] * A[i][i]
  ///   for k = i+1...n-1                 // Modify sub-matrix
  ///     for j = i+1...n-1
  ///       A[j][k] = A[j][k] – A[j][i] * A[i][k]
  ///   A[i][i] = 1 / A[i][i];            // Inverse diagonal again so that L*U = A
  /// 
  /// The pseudo-code in C++ for the corresponding sparse matrix algorithm for matrix A
  /// and separate lower (upper) triangular matrix L(U) would be:
  /// 
  ///  for 1 = 0...n-1
  ///    U[i][i] = 1 / A[i][i]
  ///    for j = i+1...n-1
  ///      L[j][i] = A[j][i] * U[i][i]
  ///    for k = i+1...n-1
  ///      for j = i+1...k
  ///        U[j][k] = A[j][k] - L[j][i] * U[i][k]
  ///      for j = k+1...n-1
  ///        L[j][k] = A[j][k] - L[j][i] * U[i][k]
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
    /// Index in U.data_ and A.data_ for U[i][i] and A[i][i] and the number of elements in the
    /// middle (j) loop for each iteration of the outer (i) loop
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> uii_aii_nj_;
    /// Number of elements in the inner (j) loops for each iteration of the middle (k) loop for the
    /// upper and lower triangular matrices
    std::vector<std::pair<std::size_t, std::size_t>> nujk_nljk_;
    /// Index in U.data_ for U[j][k], in A.data_ for A[j][k], in L.data_ for L[j][i], and
    /// in U.data_ for U[i][k] for each iteration of the inner (j) loop for the upper triangular matrix
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> ujk_ajk_lji_uik_;
    /// Index in L.data_ for L[j][k], in A.data_ for A[j][k], and in L.data_ for L[j][i], and
    /// in U.data_ for U[i][k] for each iteration of the inner (j) loop for the lower triangular matrix
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> ljk_ajk_lji_uik_;

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
    requires(SparseMatrixConcept<SparseMatrixPolicy>) LuDecompositionMozart(const SparseMatrixPolicy& matrix);

    ~LuDecompositionMozart() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>) static LuDecompositionMozart Create(const SparseMatrixPolicy& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return L and U Sparse matrices
    template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>) static std::pair<SparseMatrixPolicy, SparseMatrixPolicy> GetLUMatrices(
        const SparseMatrixPolicy& A,
        typename SparseMatrixPolicy::value_type initial_value);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param L The lower triangular matrix created by decomposition
    /// @param U The upper triangular matrix created by decomposition
    template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>) void Decompose(
        const SparseMatrixPolicy& A,
        SparseMatrixPolicy& L,
        SparseMatrixPolicy& U) const;
    template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>) void Decompose(
        const SparseMatrixPolicy& A,
        SparseMatrixPolicy& L,
        SparseMatrixPolicy& U) const;

   private:
    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>) void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_mozart.inl"
