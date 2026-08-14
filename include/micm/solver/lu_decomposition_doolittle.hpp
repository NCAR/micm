// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/types.hpp>

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
  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  class LuDecompositionDoolittle
  {
    using SparseMatrix = SparseMatrixPolicy;
    template<class U>
    using Vector = typename SparseMatrix::VectorType<U>;
    template<class U>
    using VectorView = typename SparseMatrix::VectorType<U>::ConstViewType;
   public:
    struct Views
    {
      VectorView<IndexPair> niLU_;
      VectorView<Bool> do_aik_;
      VectorView<Index> aik_;
      VectorView<IndexPair> uik_nkj_;
      VectorView<IndexPair> lij_ujk_;
      VectorView<Bool> do_aki_;
      VectorView<Index> aki_;
      VectorView<IndexPair> lki_nkj_;
      VectorView<IndexPair> lkj_uji_;
      VectorView<Index> uii_;

      Views() = default;

      Views(
        const Vector<IndexPair>& niLU,
        const Vector<Bool>& do_aik,
        const Vector<Index>& aik,
        const Vector<IndexPair>& uik_nkj,
        const Vector<IndexPair>& lij_ujk,
        const Vector<Bool>& do_aki,
        const Vector<Index>& aki,
        const Vector<IndexPair>& lki_nkj,
        const Vector<IndexPair>& lkj_uji,
        const Vector<Index>& uii
      ) : niLU_(niLU.GetView()),
          do_aik_(do_aik.GetView()),
          aik_(aik.GetView()),
          uik_nkj_(uik_nkj.GetView()),
          lij_ujk_(lij_ujk.GetView()),
          do_aki_(do_aki.GetView()),
          aki_(aki.GetView()),
          lki_nkj_(lki_nkj.GetView()),
          lkj_uji_(lkj_uji.GetView()),
          uii_(uii.GetView())
      {
      }
    };
   protected:
    /// number of elements in the middle (k) loops for lower and upper triangular matrices, respectively,
    /// for each iteration of the outer (i) loop
    Vector<IndexPair> niLU_;
    /// True when A[i][k] is non-zero for each iteration of the middle (k) loop for the upper
    /// triangular matrix; False otherwise.
    Vector<Bool> do_aik_;
    /// Index in A.data_ for A[i][k] for each iteration of the middle (k) loop for the upper
    /// triangular matrix when A[i][k] is non-zero
    Vector<Index> aik_;
    /// Index in U.data_ for U[i][k] for each iteration of the middle (k) loop for the upper
    /// triangular matrix when U[i][k] is non-zero, and the corresponding number of elements
    /// in the inner (j) loop
    Vector<IndexPair> uik_nkj_;
    /// Index in L.data_ for L[i][j], and in U.data_ for U[j][k] in the upper inner (j) loop
    /// when L[i][j] and U[j][k] are both non-zero.
    Vector<IndexPair> lij_ujk_;
    /// True when A[k][i] is non-zero for each iteration of the middle (k) loop for the lower
    /// triangular matrix; False otherwise.
    Vector<Bool> do_aki_;
    /// Index in A.data_ for A[k][i] for each iteration of the middle (k) loop for the lower
    /// triangular matrix when A[k][i] is non-zero.
    Vector<Index> aki_;
    /// Index in L.data_ for L[k][i] for each iteration of the middle (k) loop for the lower
    /// triangular matrix when L[k][i] is non-zero, and the corresponding number of elements
    /// in the inner (j) loop
    Vector<IndexPair> lki_nkj_;
    /// Index in L.data_ for L[k][j], and in U.data_ for U[j][i] in the lower inner (j) loop
    /// when L[k][j] and U[j][i] are both non-zero.
    Vector<IndexPair> lkj_uji_;
    /// Index in U.data_ for U[i][i] for each interation in the middle (k) loop for the lower
    /// triangular matrix when L[k][i] is non-zero
    Vector<Index> uii_;
    /// MICM_LAMBDA compatible views for index vectors
    Views views_;

   public:
    /// @brief default constructor
    LuDecompositionDoolittle();

    LuDecompositionDoolittle(const LuDecompositionDoolittle&) = delete;
    LuDecompositionDoolittle& operator=(const LuDecompositionDoolittle&) = delete;

    LuDecompositionDoolittle(LuDecompositionDoolittle&& other) noexcept;
    LuDecompositionDoolittle& operator=(LuDecompositionDoolittle&&) noexcept;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    LuDecompositionDoolittle(const SparseMatrixPolicy& matrix);

    ~LuDecompositionDoolittle() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    static LuDecompositionDoolittle Create(const SparseMatrixPolicy& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return L and U Sparse matrices
    static std::pair<SparseMatrixPolicy, SparseMatrixPolicy> GetLUMatrices(
        const SparseMatrixPolicy& A,
        typename SparseMatrixPolicy::value_type initial_value,
        bool indexing_only = false);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param L The lower triangular matrix created by decomposition
    /// @param U The upper triangular matrix created by decomposition
    void Decompose(const SparseMatrixPolicy& A, SparseMatrixPolicy& L, SparseMatrixPolicy& U) const;

   protected:
    /// @brief Sparse LU fill pattern of A together with the row/column adjacency
    /// needed to build the decomposition index arrays without scanning the dense
    /// (i, k, j) grid. Shared by GetLUMatrices (which only needs the id sets) and
    /// Initialize (which also walks the adjacency).
    struct FillPattern
    {
      /// Sorted non-zero positions of the L and U factors (used to build the matrices)
      std::set<std::pair<Index, Index>> L_ids_, U_ids_;
      /// Non-zero structure of the input matrix A: Arow_[r] = sorted columns,
      /// Acol_[c] = sorted rows
      std::vector<std::vector<Index>> Arow_, Acol_;
      /// Lrow_[i] = sorted columns j < i where L[i][j] != 0
      std::vector<std::vector<Index>> Lrow_;
      /// Urow_[i] = sorted columns k >= i where U[i][k] != 0
      std::vector<std::vector<Index>> Urow_;
      /// Lcol_[i] = sorted rows k > i where L[k][i] != 0
      std::vector<std::vector<Index>> Lcol_;
    };

    /// @brief Compute the sparse LU fill pattern of A in time proportional to the
    /// number of non-zeros in the factors (plus an O(n^2) scan of the sparse input A),
    /// rather than the O(n^3) dense triple loop.
    /// @param A Sparse matrix that will be decomposed
    /// @return Fill pattern and adjacency of A, L and U
    static FillPattern ComputeFillPattern(const SparseMatrixPolicy& A);

    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_doolittle.inl"
