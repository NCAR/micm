// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionMozart::LuDecompositionMozart()
  {
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart::LuDecompositionMozart(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart LuDecompositionMozart::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozart lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    MICM_PROFILE_FUNCTION();

    std::size_t n = matrix.NumRows();
    auto LU = GetLUMatrices<SparseMatrixPolicy>(matrix, initial_value);
    const auto& L_row_start = LU.first.RowStartVector();
    const auto& L_row_ids = LU.first.RowIdsVector();
    const auto& U_row_start = LU.second.RowStartVector();
    const auto& U_row_ids = LU.second.RowIdsVector();
    for (std::size_t i = 0; i < n; ++i)
    {
      std::tuple<std::size_t, std::size_t, std::size_t> lii_nuij_nlij(0, 0, 0);
      std::get<0>(lii_nuij_nlij) = LU.first.VectorIndex(0, i, i);
      // set initial values for U matrix
      for (std::size_t j = i; j < n; ++j)
      {
        if (matrix.IsZero(i, j))
        {
          if (!LU.second.IsZero(i, j))
            fill_uij_.push_back(LU.second.VectorIndex(0, i, j));
          continue;
        }
        uij_aij_.push_back(std::make_pair(LU.second.VectorIndex(0, i, j), matrix.VectorIndex(0, i, j)));
        ++(std::get<1>(lii_nuij_nlij));
      }
      // set initial values for L matrix
      for (std::size_t j = 0; j < i; ++j)
      {
        if (matrix.IsZero(i, j))
        {
          if (!LU.first.IsZero(i, j))
            fill_lij_.push_back(LU.first.VectorIndex(0, i, j));
          continue;
        }
        lij_aij_.push_back(std::make_pair(LU.first.VectorIndex(0, i, j), matrix.VectorIndex(0, i, j)));
        ++(std::get<2>(lii_nuij_nlij));
      }
      lii_nuij_nlij_.push_back(lii_nuij_nlij);
    }
    for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    {
      std::tuple<std::size_t, std::size_t, std::size_t> uii_nj_nk(0, 0, 0);
      std::get<0>(uii_nj_nk) = LU.second.VectorIndex(0, i, i);
      // middle j loop to set L[j][i]
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (LU.first.IsZero(j, i))
          continue;
        lji_.push_back(LU.first.VectorIndex(0, j, i));
        ++(std::get<1>(uii_nj_nk));
      }
      // middle k loop to set U[j][k] and L[j][k]
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (LU.second.IsZero(i, k))
          continue;
        std::tuple<std::size_t, std::size_t, std::size_t> nujk_nljk_uik(0, 0, 0);
        std::get<2>(nujk_nljk_uik) = LU.second.VectorIndex(0, i, k);
        // inner j loop to set U[j][k]
        for (std::size_t j = i + 1; j <= k; ++j)
        {
          if (LU.first.IsZero(j, i))
            continue;
          std::pair<std::size_t, std::size_t> ujk_lji(0, 0);
          ujk_lji.first = LU.second.VectorIndex(0, j, k);
          ujk_lji.second = LU.first.VectorIndex(0, j, i);
          ujk_lji_.push_back(ujk_lji);
          ++(std::get<0>(nujk_nljk_uik));
        }
        // inner j loop to set L[j][k]
        for (std::size_t j = k + 1; j < n; ++j)
        {
          if (LU.first.IsZero(j, i))
            continue;
          std::pair<std::size_t, std::size_t> ljk_lji(0, 0);
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

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LuDecompositionMozart::GetLUMatrices(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value)
  {
    MICM_PROFILE_FUNCTION();

    std::size_t n = A.NumRows();
    std::set<std::pair<std::size_t, std::size_t>> L_ids, U_ids;
    const auto& row_start = A.RowStartVector();
    const auto& row_ids = A.RowIdsVector();
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
          continue;
        for (std::size_t j = i + 1; j <= k; ++j)
          if (std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
            U_ids.insert(std::make_pair(j, k));
        for (std::size_t j = k + 1; j < n; ++j)
          if (std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
            L_ids.insert(std::make_pair(j, k));
      }
    }
    auto L_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : L_ids)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : U_ids)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LU(L_builder, U_builder);
    return LU;
  }

  template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozart::Decompose(const SparseMatrixPolicy& A, SparseMatrixPolicy& L, SparseMatrixPolicy& U)
      const
  {
    MICM_PROFILE_FUNCTION();
    const std::size_t n = A.NumRows();

    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto L_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto U_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto uij_aij = uij_aij_.begin();
      auto lij_aij = lij_aij_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();

      for (auto& lii_nuij_nlij : lii_nuij_nlij_)
      {
        for (std::size_t i = 0; i < std::get<1>(lii_nuij_nlij); ++i)
        {
          U_vector[uij_aij->first] = A_vector[uij_aij->second];
          ++uij_aij;
        }
        L_vector[std::get<0>(lii_nuij_nlij)] = 1.0;
        for (std::size_t i = 0; i < std::get<2>(lii_nuij_nlij); ++i)
        {
          L_vector[lij_aij->first] = A_vector[lij_aij->second];
          ++lij_aij;
        }
      }
      for (auto& fill_uij : fill_uij_)
        U_vector[fill_uij] = 0;
      for (auto& fill_lij : fill_lij_)
        L_vector[fill_lij] = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        auto Uii_inverse = 1.0 / U_vector[std::get<0>(*uii_nj_nk)];
        for (std::size_t ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          L_vector[*lji] = L_vector[*lji] * Uii_inverse;
          ++lji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const std::size_t uik = std::get<2>(*nujk_nljk_uik);
          for (std::size_t ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            U_vector[ujk_lji->first] -= L_vector[ujk_lji->second] * U_vector[uik];
            ++ujk_lji;
          }
          for (std::size_t ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
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
  inline void LuDecompositionMozart::Decompose(const SparseMatrixPolicy& A, SparseMatrixPolicy& L, SparseMatrixPolicy& U)
      const
  {
    MICM_PROFILE_FUNCTION();

    const std::size_t n = A.NumRows();
    const std::size_t A_BlockSize = A.NumberOfBlocks();
    constexpr std::size_t A_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t A_GroupSizeOfFlatBlockSize = A.GroupSize(A.FlatBlockSize());
    const std::size_t L_GroupSizeOfFlatBlockSize = L.GroupSize(L.FlatBlockSize());
    const std::size_t U_GroupSizeOfFlatBlockSize = U.GroupSize(U.FlatBlockSize());
    double Uii_inverse[A_GroupVectorSize];

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A_BlockSize); ++i_group)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_group * A_GroupSizeOfFlatBlockSize);
      auto L_vector = std::next(L.AsVector().begin(), i_group * L_GroupSizeOfFlatBlockSize);
      auto U_vector = std::next(U.AsVector().begin(), i_group * U_GroupSizeOfFlatBlockSize);
      auto uij_aij = uij_aij_.begin();
      auto lij_aij = lij_aij_.begin();
      auto uii_nj_nk = uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();
      const std::size_t n_cells = std::min(A_GroupVectorSize, A_BlockSize - i_group * A_GroupVectorSize);
      for (auto& lii_nuij_nlij : lii_nuij_nlij_)
      {
        for (std::size_t i = 0; i < std::get<1>(lii_nuij_nlij); ++i)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            U_vector[uij_aij->first + i_cell] = A_vector[uij_aij->second + i_cell];
          ++uij_aij;
        }
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          L_vector[std::get<0>(lii_nuij_nlij) + i_cell] = 1.0;
        for (std::size_t i = 0; i < std::get<2>(lii_nuij_nlij); ++i)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[lij_aij->first + i_cell] = A_vector[lij_aij->second + i_cell];
          ++lij_aij;
        }
      }
      for (auto& fill_uij : fill_uij_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          U_vector[fill_uij + i_cell] = 0;
      for (auto& fill_lij : fill_lij_)
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          L_vector[fill_lij + i_cell] = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          Uii_inverse[i_cell] = 1.0 / U_vector[std::get<0>(*uii_nj_nk) + i_cell];
        for (std::size_t ij = 0; ij < std::get<1>(*uii_nj_nk); ++ij)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[*lji + i_cell] = L_vector[*lji + i_cell] * Uii_inverse[i_cell];
          ++lji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(*uii_nj_nk); ++ik)
        {
          const std::size_t uik = std::get<2>(*nujk_nljk_uik);
          for (std::size_t ij = 0; ij < std::get<0>(*nujk_nljk_uik); ++ij)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[ujk_lji->first + i_cell] -= L_vector[ujk_lji->second + i_cell] * U_vector[uik + i_cell];
            ++ujk_lji;
          }
          for (std::size_t ij = 0; ij < std::get<1>(*nujk_nljk_uik); ++ij)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[ljk_lji->first + i_cell] -= L_vector[ljk_lji->second + i_cell] * U_vector[uik + i_cell];
            ++ljk_lji;
          }
          ++nujk_nljk_uik;
        }
        ++uii_nj_nk;
      }
    }
  }

}  // namespace micm
