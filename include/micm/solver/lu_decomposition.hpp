// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  class LuDecomposition
  {
    std::vector<std::pair<std::size_t, std::size_t>> niLU_;
    std::vector<bool> do_aik_;
    std::vector<std::size_t> aik_;
    std::vector<std::pair<std::size_t, std::size_t>> uik_nkj_;
    std::vector<std::pair<std::size_t, std::size_t>> lij_ujk_;
    std::vector<bool> do_aki_;
    std::vector<std::size_t> aki_;
    std::vector<std::pair<std::size_t, std::size_t>> lki_nkj_;
    std::vector<std::pair<std::size_t, std::size_t>> lkj_uji_;
    std::vector<std::size_t> uii_;

   public:
    /// @brief default constructor
    LuDecomposition();

    /// @brief Construct an LU decomposition algorithm for a given sparse matrix
    /// @param matrix Sparse matrix
    template<class T>
    LuDecomposition(const SparseMatrix<T>& matrix);

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix the will be decomposed
    /// @return L and U Sparse matrices
    template<class T>
    static std::pair<SparseMatrix<T>, SparseMatrix<T>> GetLUMatrices(const SparseMatrix<T>& A);

    /// @brief Perform an LU decomposition on a given A matrix
    /// @param A Sparse matrix to decompose
    /// @param LU the lower and upper triangular matrices returned as a square matrix
    ///           The diagonal of LU belongs to the upper triangular matrix and the
    ///           diagonal of the lower triangular matrix shoud be assumed to be 1
    template<class T>
    void Decompose(const SparseMatrix<T>& A, SparseMatrix<T>& L, SparseMatrix<T>& U) const;

    void Print()
    {
      std::cout << "niLU ";
      for (auto& elem : niLU_)
        std::cout << "( " << elem.first << ", " << elem.second << ") ";
      std::cout << std::endl;
      std::cout << "do_aik ";
      for (auto elem : do_aik_)
        std::cout << elem << " ";
      std::cout << std::endl;
      std::cout << "aik ";
      for (auto& elem : aik_)
        std::cout << elem << " ";
      std::cout << std::endl;
      std::cout << "uik_nkj ";
      for (auto& elem : uik_nkj_)
        std::cout << "( " << elem.first << ", " << elem.second << ") ";
      std::cout << std::endl;
      std::cout << "lij_ujk ";
      for (auto& elem : lij_ujk_)
        std::cout << "( " << elem.first << ", " << elem.second << ") ";
      std::cout << std::endl;
      std::cout << "do_aki ";
      for (auto elem : do_aki_)
        std::cout << elem << " ";
      std::cout << std::endl;
      std::cout << "aki ";
      for (auto& elem : aki_)
        std::cout << elem << " ";
      std::cout << std::endl;
      std::cout << "lki_nkj ";
      for (auto& elem : lki_nkj_)
        std::cout << "( " << elem.first << ", " << elem.second << ") ";
      std::cout << std::endl;
      std::cout << "lkj_uji ";
      for (auto& elem : lkj_uji_)
        std::cout << "( " << elem.first << ", " << elem.second << ") ";
      std::cout << std::endl;
      std::cout << "uii ";
      for (auto& elem : uii_)
        std::cout << elem << " ";
      std::cout << std::endl;
    }
  };

  inline LuDecomposition::LuDecomposition()
  {
  }

  template<class T>
  inline LuDecomposition::LuDecomposition(const SparseMatrix<T>& matrix)
  {
    std::size_t n = matrix[0].size();
    auto LU = GetLUMatrices(matrix);
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
          lij_ujk_.push_back(std::make_pair(j_id, LU.second.VectorIndex(0, j, k)));
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
          lkj_uji_.push_back(std::make_pair(j_id, LU.second.VectorIndex(0, j, i)));
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

  template<class T>
  inline std::pair<SparseMatrix<T>, SparseMatrix<T>> LuDecomposition::GetLUMatrices(const SparseMatrix<T>& A)
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
    auto L_builder = micm::SparseMatrix<T>::create(n).number_of_blocks(A.size());
    for (auto& pair : L_ids)
    {
      L_builder = L_builder.with_element(pair.first, pair.second);
    }
    auto U_builder = micm::SparseMatrix<T>::create(n).number_of_blocks(A.size());
    for (auto& pair : U_ids)
    {
      U_builder = U_builder.with_element(pair.first, pair.second);
    }
    std::pair<SparseMatrix<T>, SparseMatrix<T>> LU(L_builder, U_builder);
    return LU;
  }

  template<class T>
  inline void LuDecomposition::Decompose(const SparseMatrix<T>& A, SparseMatrix<T>& L, SparseMatrix<T>& U) const
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
}  // namespace micm
