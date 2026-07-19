// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  inline LuDecompositionMozart::LuDecompositionMozart() = default;

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart::LuDecompositionMozart(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart LuDecompositionMozart::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozart lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(
        matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    Index n = matrix.NumRows();
    auto LU = GetLUMatrices<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, initial_value, true);
    for (Index i = 0; i < n; ++i)
    {
      std::tuple<Index, Index, Index> lii_nuji_nlji(0, 0, 0);
      std::get<0>(lii_nuji_nlji) = LU.first.VectorIndex(0, i, i);
      // set initial values for U matrix
      for (Index j = 0; j <= i; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!LU.second.IsZero(j, i))
          {
            fill_uji_.push_back(LU.second.VectorIndex(0, j, i));
          }
          continue;
        }
        uji_aji_.push_back(std::make_pair(LU.second.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i)));
        ++(std::get<1>(lii_nuji_nlji));
      }
      // set initial values for L matrix
      for (Index j = i + 1; j < n; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!LU.first.IsZero(j, i))
          {
            fill_lji_.push_back(LU.first.VectorIndex(0, j, i));
          }
          continue;
        }
        lji_aji_.push_back(std::make_pair(LU.first.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i)));
        ++(std::get<2>(lii_nuji_nlji));
      }
      lii_nuji_nlji_.push_back(lii_nuji_nlji);
    }
    for (Index i = 0; i < matrix.NumRows(); ++i)
    {
      std::tuple<Index, Index, Index> uii_nj_nk(0, 0, 0);
      std::get<0>(uii_nj_nk) = LU.second.VectorIndex(0, i, i);
      // middle j loop to set L[j][i]
      for (Index j = i + 1; j < n; ++j)
      {
        if (LU.first.IsZero(j, i))
        {
          continue;
        }
        lji_.push_back(LU.first.VectorIndex(0, j, i));
        ++(std::get<1>(uii_nj_nk));
      }
      // middle k loop to set U[j][k] and L[j][k]
      for (Index k = i + 1; k < n; ++k)
      {
        if (LU.second.IsZero(i, k))
        {
          continue;
        }
        std::tuple<Index, Index, Index> nujk_nljk_uik(0, 0, 0);
        std::get<2>(nujk_nljk_uik) = LU.second.VectorIndex(0, i, k);
        // inner j loop to set U[j][k]
        for (Index j = i + 1; j <= k; ++j)
        {
          if (LU.first.IsZero(j, i))
          {
            continue;
          }
          std::pair<Index, Index> ujk_lji(0, 0);
          ujk_lji.first = LU.second.VectorIndex(0, j, k);
          ujk_lji.second = LU.first.VectorIndex(0, j, i);
          ujk_lji_.push_back(ujk_lji);
          ++(std::get<0>(nujk_nljk_uik));
        }
        // inner j loop to set L[j][k]
        for (Index j = k + 1; j < n; ++j)
        {
          if (LU.first.IsZero(j, i))
          {
            continue;
          }
          std::pair<Index, Index> ljk_lji(0, 0);
          ljk_lji.first = LU.first.VectorIndex(0, j, k);
          ljk_lji.second = LU.first.VectorIndex(0, j, i);
          ljk_lji_.push_back(ljk_lji);
          ++(std::get<1>(nujk_nljk_uik));
        }
        nujk_nljk_uik_.push_back(nujk_nljk_uik);
        ++(std::get<2>(uii_nj_nk));
      }
      uii_nj_nk_.push_back(uii_nj_nk);
    }
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<LMatrixPolicy, UMatrixPolicy> LuDecompositionMozart::GetLUMatrices(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    Index n = A.NumRows();
    std::set<std::pair<Index, Index>> L_ids, U_ids;
    for (Index i = 0; i < n; ++i)
    {
      for (Index j = i; j < n; ++j)
      {
        if (!A.IsZero(i, j))
        {
          U_ids.insert(std::make_pair(i, j));
        }
      }
      L_ids.insert(std::make_pair(i, i));
      for (Index j = 0; j < i; ++j)
      {
        if (!A.IsZero(i, j))
        {
          L_ids.insert(std::make_pair(i, j));
        }
      }
    }
    for (Index i = 0; i < n; ++i)
    {
      for (Index j = i + 1; j < n; ++j)
      {
        if (!A.IsZero(j, i))
        {
          L_ids.insert(std::make_pair(j, i));
        }
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (!(U_ids.contains(std::make_pair(i, k))))
        {
          continue;
        }
        for (Index j = i + 1; j <= k; ++j)
        {
          if (L_ids.contains(std::make_pair(j, i)))
          {
            U_ids.insert(std::make_pair(j, k));
          }
        }
        for (Index j = k + 1; j < n; ++j)
        {
          if (L_ids.contains(std::make_pair(j, i)))
          {
            L_ids.insert(std::make_pair(j, k));
          }
        }
      }
    }
    auto L_builder = LMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : L_ids)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = UMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : U_ids)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<LMatrixPolicy, UMatrixPolicy> LU(
        LMatrixPolicy(L_builder, indexing_only), UMatrixPolicy(U_builder, indexing_only));
    return LU;
  }

  template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    const Index n = A.NumRows();

    // Loop over blocks
    for (Index i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto L_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto U_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto uji_aji = uji_aji_.begin();
      auto lji_aji = lji_aji_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();

      for (const auto& lii_nuji_nlji : lii_nuji_nlji_)
      {
        for (Index i = 0; i < std::get<1>(lii_nuji_nlji); ++i)
        {
          U_vector[uji_aji->first] = A_vector[uji_aji->second];
          ++uji_aji;
        }
        L_vector[std::get<0>(lii_nuji_nlji)] = 1.0;
        for (Index i = 0; i < std::get<2>(lii_nuji_nlji); ++i)
        {
          L_vector[lji_aji->first] = A_vector[lji_aji->second];
          ++lji_aji;
        }
      }
      for (const auto& fill_uji : fill_uji_)
      {
        U_vector[fill_uji] = 0;
      }
      for (const auto& fill_lji : fill_lji_)
      {
        L_vector[fill_lji] = 0;
      }
      for (Index i = 0; i < n; ++i)
      {
        auto Uii_inverse = 1.0 / U_vector[std::get<0>(*uii_nj_nk)];
        for (Index ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          L_vector[*lji] = L_vector[*lji] * Uii_inverse;
          ++lji;
        }
        for (Index ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const Index uik = std::get<2>(*nujk_nljk_uik);
          for (Index ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            U_vector[ujk_lji->first] -= L_vector[ujk_lji->second] * U_vector[uik];
            ++ujk_lji;
          }
          for (Index ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
          {
            L_vector[ljk_lji->first] -= L_vector[ljk_lji->second] * U_vector[uik];
            ++ljk_lji;
          }
          ++nujk_nljk_uik;
        }
        ++uii_nj_nk;
      }
    }
  }

  template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    const Index n = A.NumRows();
    const Index A_BlockSize = A.NumberOfBlocks();
    constexpr Index A_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const Index A_GroupSizeOfFlatBlockSize = A.GroupSize();
    const Index L_GroupSizeOfFlatBlockSize = L.GroupSize();
    const Index U_GroupSizeOfFlatBlockSize = U.GroupSize();
    Real Uii_inverse[A_GroupVectorSize];

    // Loop over groups of blocks
    for (Index i_group = 0; i_group < A.NumberOfGroups(A_BlockSize); ++i_group)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_group * A_GroupSizeOfFlatBlockSize);
      auto L_vector = std::next(L.AsVector().begin(), i_group * L_GroupSizeOfFlatBlockSize);
      auto U_vector = std::next(U.AsVector().begin(), i_group * U_GroupSizeOfFlatBlockSize);
      auto uji_aji = uji_aji_.begin();
      auto lji_aji = lji_aji_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();
      const Index n_cells = std::min(A_GroupVectorSize, A_BlockSize - i_group * A_GroupVectorSize);
      for (const auto& lii_nuji_nlji : lii_nuji_nlji_)
      {
        for (Index i = 0; i < std::get<1>(lii_nuji_nlji); ++i)
        {
          for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            U_vector[uji_aji->first + i_cell] = A_vector[uji_aji->second + i_cell];
          }
          ++uji_aji;
        }
        for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          L_vector[std::get<0>(lii_nuji_nlji) + i_cell] = 1.0;
        }
        for (Index i = 0; i < std::get<2>(lii_nuji_nlji); ++i)
        {
          for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            L_vector[lji_aji->first + i_cell] = A_vector[lji_aji->second + i_cell];
          }
          ++lji_aji;
        }
      }
      for (const auto& fill_uji : fill_uji_)
      {
        for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          U_vector[fill_uji + i_cell] = 0;
        }
      }
      for (const auto& fill_lji : fill_lji_)
      {
        for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          L_vector[fill_lji + i_cell] = 0;
        }
      }
      for (Index i = 0; i < n; ++i)
      {
        for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          Uii_inverse[i_cell] = 1.0 / U_vector[std::get<0>(*uii_nj_nk) + i_cell];
        }
        for (Index ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            L_vector[*lji + i_cell] = L_vector[*lji + i_cell] * Uii_inverse[i_cell];
          }
          ++lji;
        }
        for (Index ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const Index uik = std::get<2>(*nujk_nljk_uik);
          for (Index ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              U_vector[ujk_lji->first + i_cell] -= L_vector[ujk_lji->second + i_cell] * U_vector[uik + i_cell];
            }
            ++ujk_lji;
          }
          for (Index ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
          {
            for (Index i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              L_vector[ljk_lji->first + i_cell] -= L_vector[ljk_lji->second + i_cell] * U_vector[uik + i_cell];
            }
            ++ljk_lji;
          }
          ++nujk_nljk_uik;
        }
        ++uii_nj_nk;
      }
    }
  }

}  // namespace micm
