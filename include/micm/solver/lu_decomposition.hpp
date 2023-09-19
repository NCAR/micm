// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix
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
  class LuDecomposition
  {
    public:
    /// number of elements in the middle (k) loops for lower and upper triangular matrices, respectively,
    /// for each iteration of the outer (i) loop
    std::vector<std::pair<std::size_t, std::size_t>> niLU_;
    /// True when A[i][k] is non-zero for each iteration of the middle (k) loop for the upper
    /// triangular matrix; False otherwise
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
    /// triangular matrix; False otherwise
    std::vector<char> do_aki_;
    /// Index in A.data_ for A[k][i] for each iteration of the middle (k) loop for the lower
    /// triangular matrix when A[k][i] is non-zero
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
    LuDecomposition();

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<typename T, typename OrderingPolicy>
    LuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix the will be decomposed
    /// @return L and U Sparse matrices
    template<typename T, typename OrderingPolicy>
    static std::pair<SparseMatrix<T, OrderingPolicy>, SparseMatrix<T, OrderingPolicy>> GetLUMatrices(
        const SparseMatrix<T, OrderingPolicy>& A,
        T initial_value);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param LU the lower and upper triangular matrices returned as a square matrix
    ///           The diagonal of LU belongs to the upper triangular matrix and the
    ///           diagonal of the lower triangular matrix shoud be assumed to be 1
    template<typename T, template<class> class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy<T>>) void Decompose(
        const SparseMatrixPolicy<T>& A,
        SparseMatrixPolicy<T>& L,
        SparseMatrixPolicy<T>& U) const;
    template<typename T, template<class> class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy<T>>) void Decompose(
        const SparseMatrixPolicy<T>& A,
        SparseMatrixPolicy<T>& L,
        SparseMatrixPolicy<T>& U) const;
  };

  inline LuDecomposition::LuDecomposition()
  {
  }

  template<typename T, typename OrderingPolicy>
  inline LuDecomposition::LuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix)
  {
    std::size_t n = matrix[0].size();
    auto LU = GetLUMatrices(matrix, T{});
    const auto& L_row_start = LU.first.RowStartVector();
    const auto& L_row_ids = LU.first.RowIdsVector();
    const auto& U_row_start = LU.second.RowStartVector();
    const auto& U_row_ids = LU.second.RowIdsVector();
    for (std::size_t i = 0; i < matrix[0].size(); ++i)
    {
      std::pair<std::size_t, std::size_t> iLU(0, 0);
      // Upper triangular matrix
      for (std::size_t k = i; k < n; ++k)
      {
        std::size_t nkj = 0;
        for (std::size_t j_id = L_row_start[i]; j_id < L_row_start[i + 1]; ++j_id)
        {
          std::size_t j = L_row_ids[j_id];
          if (j >= i)
            break;
          if (LU.second.IsZero(j, k))
            continue;
          ++nkj;
          lij_ujk_.push_back(std::make_pair(LU.first.VectorIndex(0, i, j), LU.second.VectorIndex(0, j, k)));
        }
        if (matrix.IsZero(i, k))
        {
          if (nkj == 0 && k != i)
            continue;
          do_aik_.push_back(false);
        }
        else
        {
          do_aik_.push_back(true);
          aik_.push_back(matrix.VectorIndex(0, i, k));
        }
        uik_nkj_.push_back(std::make_pair(LU.second.VectorIndex(0, i, k), nkj));
        ++(iLU.second);
      }
      // Lower triangular matrix
      lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, i, i), 0));
      for (std::size_t k = i + 1; k < n; ++k)
      {
        std::size_t nkj = 0;
        for (std::size_t j_id = L_row_start[k]; j_id < L_row_start[k + 1]; ++j_id)
        {
          std::size_t j = L_row_ids[j_id];
          if (j >= i)
            break;
          if (LU.second.IsZero(j, i))
            continue;
          ++nkj;
          lkj_uji_.push_back(std::make_pair(LU.first.VectorIndex(0, k, j), LU.second.VectorIndex(0, j, i)));
        }
        if (matrix.IsZero(k, i))
        {
          if (nkj == 0)
            continue;
          do_aki_.push_back(false);
        }
        else
        {
          do_aki_.push_back(true);
          aki_.push_back(matrix.VectorIndex(0, k, i));
        }
        uii_.push_back(LU.second.VectorIndex(0, i, i));
        lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, k, i), nkj));
        ++(iLU.first);
      }
      niLU_.push_back(iLU);
    }
  }

  template<typename T, typename OrderingPolicy>
  inline std::pair<SparseMatrix<T, OrderingPolicy>, SparseMatrix<T, OrderingPolicy>> LuDecomposition::GetLUMatrices(
      const SparseMatrix<T, OrderingPolicy>& A,
      T initial_value)
  {
    std::size_t n = A[0].size();
    std::set<std::pair<std::size_t, std::size_t>> L_ids, U_ids;
    const auto& row_start = A.RowStartVector();
    const auto& row_ids = A.RowIdsVector();
    for (std::size_t i = 0; i < n; ++i)
    {
      // Upper triangular matrix
      for (std::size_t k = i; k < n; ++k)
      {
        if (!A.IsZero(i, k) || k == i)
        {
          U_ids.insert(std::make_pair(i, k));
          continue;
        }
        for (std::size_t j = 0; j < i; ++j)
        {
          if (L_ids.find(std::make_pair(i, j)) != L_ids.end() && U_ids.find(std::make_pair(j, k)) != U_ids.end())
          {
            U_ids.insert(std::make_pair(i, k));
            break;
          }
        }
      }
      // Lower triangular matrix
      for (std::size_t k = i; k < n; ++k)
      {
        if (!A.IsZero(k, i) || k == i)
        {
          L_ids.insert(std::make_pair(k, i));
          continue;
        }
        for (std::size_t j = 0; j < i; ++j)
        {
          if (L_ids.find(std::make_pair(k, j)) != L_ids.end() && U_ids.find(std::make_pair(j, i)) != U_ids.end())
          {
            L_ids.insert(std::make_pair(k, i));
            break;
          }
        }
      }
    }
    auto L_builder =
        micm::SparseMatrix<T, OrderingPolicy>::create(n).number_of_blocks(A.size()).initial_value(initial_value);
    for (auto& pair : L_ids)
    {
      L_builder = L_builder.with_element(pair.first, pair.second);
    }
    auto U_builder =
        micm::SparseMatrix<T, OrderingPolicy>::create(n).number_of_blocks(A.size()).initial_value(initial_value);
    for (auto& pair : U_ids)
    {
      U_builder = U_builder.with_element(pair.first, pair.second);
    }
    std::pair<SparseMatrix<T, OrderingPolicy>, SparseMatrix<T, OrderingPolicy>> LU(L_builder, U_builder);
    return LU;
  }

  template<typename T, template<class> class SparseMatrixPolicy>
  requires(!VectorizableSparse<SparseMatrixPolicy<T>>) void LuDecomposition::Decompose(
      const SparseMatrixPolicy<T>& A,
      SparseMatrixPolicy<T>& L,
      SparseMatrixPolicy<T>& U) const
  {
    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.size(); ++i_block)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto L_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto U_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto do_aik = do_aik_.begin();
      auto aik = aik_.begin();
      auto uik_nkj = uik_nkj_.begin();
      auto lij_ujk = lij_ujk_.begin();
      auto do_aki = do_aki_.begin();
      auto aki = aki_.begin();
      auto lki_nkj = lki_nkj_.begin();
      auto lkj_uji = lkj_uji_.begin();
      auto uii = uii_.begin();
      for (auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          if (*(do_aik++))
            U_vector[uik_nkj->first] = A_vector[*(aik++)];
          for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
          {
            U_vector[uik_nkj->first] -= L_vector[lij_ujk->first] * U_vector[lij_ujk->second];
            ++lij_ujk;
          }
          ++uik_nkj;
        }
        // Lower triangular matrix
        L_vector[(lki_nkj++)->first] = 1.0;
        for (std::size_t iL = 0; iL < inLU.first; ++iL)
        {
          if (*(do_aki++))
            L_vector[lki_nkj->first] = A_vector[*(aki++)];
          for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
          {
            L_vector[lki_nkj->first] -= L_vector[lkj_uji->first] * U_vector[lkj_uji->second];
            ++lkj_uji;
          }
          L_vector[lki_nkj->first] /= U_vector[*uii];
          ++lki_nkj;
          ++uii;
        }
      }
    }
  }

  template<typename T, template<class> class SparseMatrixPolicy>
  requires(VectorizableSparse<SparseMatrixPolicy<T>>) void LuDecomposition::Decompose(
      const SparseMatrixPolicy<T>& A,
      SparseMatrixPolicy<T>& L,
      SparseMatrixPolicy<T>& U) const
  {
    const std::size_t n_cells = A.GroupVectorSize();
    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A.size()); ++i_group)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_group * A.GroupSize(A.FlatBlockSize()));
      auto L_vector = std::next(L.AsVector().begin(), i_group * L.GroupSize(L.FlatBlockSize()));
      auto U_vector = std::next(U.AsVector().begin(), i_group * U.GroupSize(U.FlatBlockSize()));
      auto do_aik = do_aik_.begin();
      auto aik = aik_.begin();
      auto uik_nkj = uik_nkj_.begin();
      auto lij_ujk = lij_ujk_.begin();
      auto do_aki = do_aki_.begin();
      auto aki = aki_.begin();
      auto lki_nkj = lki_nkj_.begin();
      auto lkj_uji = lkj_uji_.begin();
      auto uii = uii_.begin();
      for (auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          if (*(do_aik++))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell){
              int U_index = uik_nkj->first + i_cell; 
              int A_index = *aik + i_cell;
              U_vector[U_index] = A_vector[A_index];
              std::cout << "this is cpu U_index: "<<U_index<<std::endl; 
              std::cout << "This is cpu A_index:"<<A_index<<std::endl; 
              std::cout << "this is cpu U value:"<< U_vector[U_index]<<std::endl; }
            ++aik;
          }
          for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[uik_nkj->first + i_cell] -= L_vector[lij_ujk->first + i_cell] * U_vector[lij_ujk->second + i_cell];
            ++lij_ujk;
          }
          ++uik_nkj;
        }
        // Lower triangular matrix
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          L_vector[lki_nkj->first + i_cell] = 1.0;
        ++lki_nkj;
        for (std::size_t iL = 0; iL < inLU.first; ++iL)
        {
          if (*(do_aki++))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[lki_nkj->first + i_cell] = A_vector[*aki + i_cell];
            ++aki;
          }
          for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[lki_nkj->first + i_cell] -= L_vector[lkj_uji->first + i_cell] * U_vector[lkj_uji->second + i_cell];
            ++lkj_uji;
          }
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[lki_nkj->first + i_cell] /= U_vector[*uii + i_cell];
          ++lki_nkj;
          ++uii;
        }
      }
    }
  }
}  // namespace micm
