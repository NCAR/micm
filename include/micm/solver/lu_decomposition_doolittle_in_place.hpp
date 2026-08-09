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
  /// The sudo-code for the corresponding dense matrix algorithm for matrix A (in-line) would be:
  ///
  /// for i = 0...n-1                 // Outer loop over rows (columns) for upper (lower) triangular matrix
  ///   for k = i...n-1               // Middle loop over columns for upper triangular matrix
  ///     for j = 0...i-1             // Inner loop over columns (rows) for lower (upper) triangular matrix
  ///       A[i][k] -= A[i][j] * A[j][k]
  ///   for k = i+1...n-1             // Middle loop over rows for lower triangular matrix
  ///     for j = 0...i-1             // Inner loop over columns (rows) for lower (upper) triangular matrix
  ///       A[k][i] -= A[k][j] * A[j][i];
  ///     A[k][i] /= A[i][i]
  ///
  /// For the sparse matrix algorithm, the indices of non-zero terms are stored in
  /// several arrays during construction. These arrays are iterated through during
  /// calls to Decompose to do the actual decomposition.
  /// Our LU Decomposition only assigns the values of the jacobian to the LU matrices
  /// when the *jacobian* is nonzero. However, the sparsity pattern of the jacobian doesn't
  /// necessarily match that of the LU matrices. There can be more nonzero elements in the LU matrices
  /// than in the jacobian. It is expected that the elements of the L and U matrices that are zero in the A matrix
  /// will be set to zero before the combined matrix is passed to the decomposition function.
  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  class LuDecompositionDoolittleInPlace
  {
    using SparseMatrix = SparseMatrixPolicy;
    template<class U>
    using Vector = typename SparseMatrix::VectorType<U>;
    template<class U>
    using VectorView = typename SparseMatrix::VectorType<U>::ConstViewType;
   public:
    struct IndexTrio
    {
      Index first_;
      Index second_;
      Index third_;
    };
    struct IndexPair
    {
      Index first_;
      Index second_;
    };
    struct Views
    {
      VectorView<IndexTrio> nik_nki_aii_;
      VectorView<IndexPair> aik_njk_;
      VectorView<IndexPair> aij_ajk_;
      VectorView<IndexPair> aki_nji_;
      VectorView<IndexPair> akj_aji_;

      Views() = default;

      Views(
        const Vector<IndexTrio>& nik_nki_aii,
        const Vector<IndexPair>& aik_njk,
        const Vector<IndexPair>& aij_ajk,
        const Vector<IndexPair>& aki_nji,
        const Vector<IndexPair>& akj_aji
      ) : nik_nki_aii_(nik_nki_aii.GetView()),
          aik_njk_(aik_njk.GetView()),
          aij_ajk_(aij_ajk.GetView()),
          aki_nji_(aki_nji.GetView()),
          akj_aji_(akj_aji.GetView())
      {
      }
    };    
   protected:
    /// Number of elements in the middle (k) loops for lower and upper triangular matrices, respectively,
    /// and the index in A.data_ for A[i][i] for each iteration of the outer (i) loop
    Vector<IndexTrio> nik_nki_aii_;
    /// Index in A.data_ for A[i][k] for each iteration of the upper middle (k) loop and the
    /// number of elements in the inner (j) loop for each upper (k) element used to set A[i][k]
    Vector<IndexPair> aik_njk_;
    /// Index in A.data_ for A[i][j] and A[j][k] for each iteration of the upper inner (j) loop
    Vector<IndexPair> aij_ajk_;
    /// Index in A.data_ for A[k][i] for each iteration of the lower middle (k) loop, the
    /// number of elements in the inner (j) loop for each lower (k) element used to set A[k][i]
    Vector<IndexPair> aki_nji_;
    /// Index in A.data_ for A[k][j] and A[j][i] for each iteration of the lower inner (j) loop
    Vector<IndexPair> akj_aji_;
    /// MICM_LAMBDA compatible views for index vectors
    Views views_;

   public:
    /// @brief default constructor
    LuDecompositionDoolittleInPlace();

    LuDecompositionDoolittleInPlace(const LuDecompositionDoolittleInPlace&) = delete;
    LuDecompositionDoolittleInPlace& operator=(const LuDecompositionDoolittleInPlace&) = delete;

    LuDecompositionDoolittleInPlace(LuDecompositionDoolittleInPlace&& other) noexcept;
    LuDecompositionDoolittleInPlace& operator=(LuDecompositionDoolittleInPlace&&) noexcept;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    LuDecompositionDoolittleInPlace(const SparseMatrixPolicy& matrix);

    ~LuDecompositionDoolittleInPlace() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    static LuDecompositionDoolittleInPlace Create(const SparseMatrixPolicy& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return L and U Sparse matrices
    static SparseMatrixPolicy GetLUMatrix(
        const SparseMatrixPolicy& A,
        typename SparseMatrixPolicy::value_type initial_value,
        bool indexing_only = false);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param L The lower triangular matrix created by decomposition
    /// @param U The upper triangular matrix created by decomposition
    void Decompose(SparseMatrixPolicy& ALU) const;

   protected:
    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_doolittle_in_place.inl"
