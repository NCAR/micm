// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/profiler/instrumentation.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix
  ///
  /// This LU decomposition uses the algorithm from the CHEM_PREPROCESSOR at: 
  /// https://github.com/ESCOMP/CHEM_PREPROCESSOR/blob/f923081508f4264e61fcef2753a9898e55d1598e/src/cam_chempp/make_lu_fac.f#L127-L214
  ///
  /// The pseudo-code in C++ for the corresponding sparse matrix algorithm for matrix A (inline change) would be: 
  ///
  /// for i = 0...n-1                       // Outer loop over columns of the sparse matrix A
  ///     if (A[i][i] != 0)
  ///     {
  ///         A[i][i] = 1 / A[i][i]         // Form diagonal inverse
  ///     }
  ///     for j = i+1...n-1                 // Multiply column below diagonal
  ///         if (A[j][i] != 0)
  ///         {
  ///             A[j][i] = A[j][i] * A[i][i]
  ///         } 
  ///     for k = i+1...n-1                 // Modify sub-matrix
  ///         if (A[i][k] != 0)
  ///         {
  ///             for j = i+1...n-1
  ///                 if (A[j][i] != 0)
  ///                 {
  ///                     if (A[j][k] != 0)
  ///                     {
  ///                         A[j][k] = A[j][k] – A[j][i] * A[i][k]
  ///                     }                     
  ///                     else
  ///                     {
  ///                         A[j][k] = - A[j][i] * A[i][k]
  ///                     }
  ///                 }
  ///         }
  ///     A[i][i] = 1 / A[i][i];            // Inverse diagonal again so that L*U = A
  /// 
  /// The pseudo-code in C++ for the corresponding sparse matrix algorithm for matrix A
  /// and lower (upper) triangular matrix L(U) would be:
  /// 
  /// for j = 1...n-1                       // Initialize the L matrix
  ///     for k = 0...j-1
  ///         L[j][k] = A[j][k]
  /// for k = 0...n-1                       // Initialize the U matrix
  ///     for j = k...n-1
  ///         U[j][k] = A[j][k]
  /// for i = 0...n-1                       // Outer loop over columns of the sparse matrix A
  ///     denom = 1 / U[i][i]               // Save the inverse of diagonal (may check if U[i][i] equals to zero or not)
  ///     for j = i+1...n-1                 // Multiply column below diagonal
  ///         if (L[j][i] != 0)
  ///         {
  ///             L[j][i] = L[j][i] * denom
  ///         } 
  ///     for k = i+1...n-1                 // Modify sub-matrix
  ///         if (U[i][k] != 0)
  ///         {
  ///             for j = i+1...k           // Modify upper triangular matrix
  ///                 if (L[j][i] != 0)
  ///                 {
  ///                     if (U[j][k] != 0)
  ///                     {
  ///                         U[j][k] = U[j][k] – L[j][i] * U[i][k]
  ///                     }                     
  ///                     else
  ///                     {
  ///                         U[j][k] = - L[j][i] * U[i][k]
  ///                     }
  ///                 }
  ///             for j = k+1...n-1         // Modify lower triangular matrix
  ///                 if (L[j][i] != 0)
  ///                 {
  ///                     if (L[j][k] != 0)
  ///                     {
  ///                         L[j][k] = L[j][k] – L[j][i] * U[i][k]
  ///                     }                     
  ///                     else
  ///                     {
  ///                         L[j][k] = - L[j][i] * U[i][k]
  ///                     }
  ///                 }
  ///         }
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
  class ChempreprocessorLuDecomposition
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
    ChempreprocessorLuDecomposition();

    ChempreprocessorLuDecomposition(const ChempreprocessorLuDecomposition&) = delete;
    ChempreprocessorLuDecomposition& operator=(const ChempreprocessorLuDecomposition&) = delete;

    ChempreprocessorLuDecomposition(ChempreprocessorLuDecomposition&& other) = default;
    ChempreprocessorLuDecomposition& operator=(ChempreprocessorLuDecomposition&&) = default;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>) ChempreprocessorLuDecomposition(const SparseMatrixPolicy& matrix);

    ~ChempreprocessorLuDecomposition() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>) static ChempreprocessorLuDecomposition
        Create(const SparseMatrixPolicy& matrix);

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

#include "chempreprocessor_lu_decomposition.inl"
