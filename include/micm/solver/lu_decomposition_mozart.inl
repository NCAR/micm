// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionMozart::LuDecompositionMozart()
  {
  }

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
    std::size_t n = matrix.NumRows();
    auto lu = GetLUMatrices<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, initial_value, true);
    for (std::size_t i = 0; i < n; ++i)
    {
      std::tuple<std::size_t, std::size_t, std::size_t> lii_nuji_nlji(0, 0, 0);
      std::get<0>(lii_nuji_nlji) = lu.first.VectorIndex(0, i, i);
      // set initial values for U matrix
      for (std::size_t j = 0; j <= i; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!lu.second.IsZero(j, i))
            fill_uji_.push_back(lu.second.VectorIndex(0, j, i));
          continue;
        }
        uji_aji_.push_back(std::make_pair(lu.second.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i)));
        ++(std::get<1>(lii_nuji_nlji));
      }
      // set initial values for L matrix
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!lu.first.IsZero(j, i))
            fill_lji_.push_back(lu.first.VectorIndex(0, j, i));
          continue;
        }
        lji_aji_.push_back(std::make_pair(lu.first.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i)));
        ++(std::get<2>(lii_nuji_nlji));
      }
      lii_nuji_nlji_.push_back(lii_nuji_nlji);
    }
    for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    {
      std::tuple<std::size_t, std::size_t, std::size_t> uii_nj_nk(0, 0, 0);
      std::get<0>(uii_nj_nk) = lu.second.VectorIndex(0, i, i);
      // middle j loop to set L[j][i]
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (lu.first.IsZero(j, i))
        {
          continue;
        }
        lji_.push_back(lu.first.VectorIndex(0, j, i));
        ++(std::get<1>(uii_nj_nk));
      }
      // middle k loop to set U[j][k] and L[j][k]
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (lu.second.IsZero(i, k))
        {
          continue;
        }
        std::tuple<std::size_t, std::size_t, std::size_t> nujk_nljk_uik(0, 0, 0);
        std::get<2>(nujk_nljk_uik) = lu.second.VectorIndex(0, i, k);
        // inner j loop to set U[j][k]
        for (std::size_t j = i + 1; j <= k; ++j)
        {
          if (lu.first.IsZero(j, i))
          {
            continue;
          }
          std::pair<std::size_t, std::size_t> ujk_lji(0, 0);
          ujk_lji.first = lu.second.VectorIndex(0, j, k);
          ujk_lji.second = lu.first.VectorIndex(0, j, i);
          ujk_lji_.push_back(ujk_lji);
          ++(std::get<0>(nujk_nljk_uik));
        }
        // inner j loop to set L[j][k]
        for (std::size_t j = k + 1; j < n; ++j)
        {
          if (lu.first.IsZero(j, i))
          {
            continue;
          }
          std::pair<std::size_t, std::size_t> ljk_lji(0, 0);
          ljk_lji.first = lu.first.VectorIndex(0, j, k);
          ljk_lji.second = lu.first.VectorIndex(0, j, i);
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
    std::size_t n = A.NumRows();
    std::set<std::pair<std::size_t, std::size_t>> L_ids, U_ids;
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = i; j < n; ++j)
        if (!A.IsZero(i, j))
          U_ids.insert(std::make_pair(i, j));
      L_ids.insert(std::make_pair(i, i));
      for (std::size_t j = 0; j < i; ++j)
        if (!A.IsZero(i, j))
          L_ids.insert(std::make_pair(i, j));
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (!A.IsZero(j, i))
          L_ids.insert(std::make_pair(j, i));
      }
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (!(std::find(U_ids.begin(), U_ids.end(), std::make_pair(i, k)) != U_ids.end()))
        {
          continue;
        }
        for (std::size_t j = i + 1; j <= k; ++j)
          if (std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
            U_ids.insert(std::make_pair(j, k));
        for (std::size_t j = k + 1; j < n; ++j)
          if (std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
            L_ids.insert(std::make_pair(j, k));
      }
    }
    auto l_builder = LMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : L_ids)
    {
      l_builder = l_builder.WithElement(pair.first, pair.second);
    }
    auto u_builder = UMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : U_ids)
    {
      u_builder = u_builder.WithElement(pair.first, pair.second);
    }
    std::pair<LMatrixPolicy, UMatrixPolicy> LU(
        LMatrixPolicy(l_builder, indexing_only), UMatrixPolicy(u_builder, indexing_only));
    return LU;
  }

  template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    const std::size_t N = A.NumRows();

    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
    {
      auto a_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto l_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto u_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto uji_aji = uji_aji_.begin();
      auto lji_aji = lji_aji_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();

      for (auto& lii_nuji_nlji : lii_nuji_nlji_)
      {
        for (std::size_t i = 0; i < std::get<1>(lii_nuji_nlji); ++i)
        {
          u_vector[uji_aji->first] = a_vector[uji_aji->second];
          ++uji_aji;
        }
        l_vector[std::get<0>(lii_nuji_nlji)] = 1.0;
        for (std::size_t i = 0; i < std::get<2>(lii_nuji_nlji); ++i)
        {
          l_vector[lji_aji->first] = a_vector[lji_aji->second];
          ++lji_aji;
        }
      }
      for (auto& fill_uji : fill_uji_)
        u_vector[fill_uji] = 0;
      for (auto& fill_lji : fill_lji_)
        l_vector[fill_lji] = 0;
      for (std::size_t i = 0; i < N; ++i)
      {
        auto uii_inverse = 1.0 / u_vector[std::get<0>(*uii_nj_nk)];
        for (std::size_t ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          l_vector[*lji] = l_vector[*lji] * uii_inverse;
          ++lji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const std::size_t UIK = std::get<2>(*nujk_nljk_uik);
          for (std::size_t ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            u_vector[ujk_lji->first] -= l_vector[ujk_lji->second] * u_vector[UIK];
            ++ujk_lji;
          }
          for (std::size_t ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
          {
            l_vector[ljk_lji->first] -= l_vector[ljk_lji->second] * u_vector[UIK];
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
    const std::size_t N = A.NumRows();
    const std::size_t A_BLOCK_SIZE = A.NumberOfBlocks();
    constexpr std::size_t A_GROUP_VECTOR_SIZE = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t A_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = A.GroupSize();
    const std::size_t L_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = L.GroupSize();
    const std::size_t U_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = U.GroupSize();
    double uii_inverse[A_GROUP_VECTOR_SIZE];

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A_BLOCK_SIZE); ++i_group)
    {
      auto a_vector = std::next(A.AsVector().begin(), i_group * A_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto l_vector = std::next(L.AsVector().begin(), i_group * L_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto u_vector = std::next(U.AsVector().begin(), i_group * U_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto uji_aji = uji_aji_.begin();
      auto lji_aji = lji_aji_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();
      const std::size_t N_CELLS = std::min(A_GROUP_VECTOR_SIZE, A_BLOCK_SIZE - i_group * A_GROUP_VECTOR_SIZE);
      for (auto& lii_nuji_nlji : lii_nuji_nlji_)
      {
        for (std::size_t i = 0; i < std::get<1>(lii_nuji_nlji); ++i)
        {
          for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
            u_vector[uji_aji->first + i_cell] = a_vector[uji_aji->second + i_cell];
          ++uji_aji;
        }
        for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          l_vector[std::get<0>(lii_nuji_nlji) + i_cell] = 1.0;
        for (std::size_t i = 0; i < std::get<2>(lii_nuji_nlji); ++i)
        {
          for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
            l_vector[lji_aji->first + i_cell] = a_vector[lji_aji->second + i_cell];
          ++lji_aji;
        }
      }
      for (auto& fill_uji : fill_uji_)
        for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          u_vector[fill_uji + i_cell] = 0;
      for (auto& fill_lji : fill_lji_)
        for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          l_vector[fill_lji + i_cell] = 0;
      for (std::size_t i = 0; i < N; ++i)
      {
        for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          uii_inverse[i_cell] = 1.0 / u_vector[std::get<0>(*uii_nj_nk) + i_cell];
        for (std::size_t ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          {
            l_vector[*lji + i_cell] = l_vector[*lji + i_cell] * uii_inverse[i_cell];
          }
          ++lji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const std::size_t UIK = std::get<2>(*nujk_nljk_uik);
          for (std::size_t ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
            {
              u_vector[ujk_lji->first + i_cell] -= l_vector[ujk_lji->second + i_cell] * u_vector[UIK + i_cell];
            }
            ++ujk_lji;
          }
          for (std::size_t ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
          {
            for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
            {
              l_vector[ljk_lji->first + i_cell] -= l_vector[ljk_lji->second + i_cell] * u_vector[UIK + i_cell];
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
