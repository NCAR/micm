// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
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
  /// For the sparse matrix algorithm, the indices of non-zero terms are stored in
  /// several arrays during construction. These arrays are iterated through during
  /// calls to Decompose to do the actual decomposition.
  ///
  /// The GetLUMatrices function creates a new sparse matrix that includes the superset
  /// of the non-zero elements in the L and U matrices. It is expected that the elements
  /// of the L and U matrices that are zero in the A matrix will be set to zero before the
  /// combined matrix is passed to the decomposition function.
  class LuDecompositionMozartInPlace
  {
   protected:
    /// Index in A.data_ for all diagonal elements, the number of iterations of the inner (j) loop
    /// for each (i) used to set A[j][i], and the number of iterations of the middle (k) loop for
    /// each (i) used to set A[j][k]
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> aii_nji_nki_;
    /// Index in A.data_ for A[j][i] for each iteration of the inner (j) loop
    /// used to set the value of A[j][i]
    std::vector<std::size_t> aji_;
    /// Index in A.data_ for A[i][k] for each iteration of the middle (k) loop,
    /// and the number of iterations of the inner (j) loop for each (k) used to set A[j][k]
    std::vector<std::pair<std::size_t, std::size_t>> aik_njk_;
    /// Index in A.data_ for A[j][k] and A[j][i] for each iteration of the inner (j) loop
    std::vector<std::pair<std::size_t, std::size_t>> ajk_aji_;

   public:
    /// @brief default constructor
    LuDecompositionMozartInPlace();

    LuDecompositionMozartInPlace(const LuDecompositionMozartInPlace&) = delete;
    LuDecompositionMozartInPlace& operator=(const LuDecompositionMozartInPlace&) = delete;

    LuDecompositionMozartInPlace(LuDecompositionMozartInPlace&& other) = default;
    LuDecompositionMozartInPlace& operator=(LuDecompositionMozartInPlace&&) = default;

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    LuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix);

    ~LuDecompositionMozartInPlace() = default;

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static LuDecompositionMozartInPlace Create(const SparseMatrixPolicy& matrix);

    /// @brief Create a combined sparse L and U matrix for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @return combined L and U Sparse matrices
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static SparseMatrixPolicy GetLUMatrix(
        const SparseMatrixPolicy& A,
        typename SparseMatrixPolicy::value_type initial_value,
        std::size_t number_of_grid_cells,
        bool empty_matrix = false);

    /// @brief Perform an LU decomposition on a given A matrix.
    ///        All elements of L and U that are zero in A should be set to zero before calling this function.
    /// @param ALU Sparse matrix to decompose (will be overwritten with L and U matrices)
    template<class SparseMatrixPolicy>
      requires(!VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(SparseMatrixPolicy& ALU) const;
    template<class SparseMatrixPolicy>
      requires(VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(SparseMatrixPolicy& ALU) const;

   protected:
    /// @brief Initialize arrays for the LU decomposition
    /// @param A Sparse matrix to decompose
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    void Initialize(const SparseMatrixPolicy& matrix, auto initial_value);
  };

}  // namespace micm

#include "lu_decomposition_mozart_in_place.inl"
