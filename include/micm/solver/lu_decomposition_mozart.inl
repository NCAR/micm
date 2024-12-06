// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionMozart::LuDecompositionMozart()
  {
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline LuDecompositionMozart::LuDecompositionMozart(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline LuDecompositionMozart LuDecompositionMozart::Create(
      const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozart lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline void LuDecompositionMozart::Initialize(
      const SparseMatrixPolicy& matrix,
      auto initial_value)
  {
    MICM_PROFILE_FUNCTION();

    std::size_t n = matrix.NumRows();
    auto LU = GetLUMatrices<SparseMatrixPolicy>(matrix, initial_value);
    const auto& L_row_start = LU.first.RowStartVector();
    const auto& L_row_ids = LU.first.RowIdsVector();
    const auto& U_row_start = LU.second.RowStartVector();
    const auto& U_row_ids = LU.second.RowIdsVector();
    l00_u00_a00_ = std::tuple(LU.first.VectorIndex(0, 0, 0),
                              LU.second.VectorIndex(0, 0, 0),
                              matrix.VectorIndex(0, 0, 0));
    for (std::size_t i = 1; i < n; ++i)
    {
      if (LU.first.IsZero(i, 0))
        continue;
      std::tuple<std::size_t, std::size_t, char> li0_ai0_doFill(0, 0, true);
      std::get<0>(li0_ai0_doFill) = LU.first.VectorIndex(0, i, 0);
      if (!matrix.IsZero(i, 0))
      {
        std::get<1>(li0_ai0_doFill) = matrix.VectorIndex(0, i, 0);
        std::get<2>(li0_ai0_doFill) = false;
      }
      li0_ai0_doFill_.push_back(li0_ai0_doFill);
    }
    for (std::size_t i = 1; i < n; ++i)
    {
      if (LU.second.IsZero(0, i))
        continue;
      std::tuple<std::size_t, std::size_t, char> u0i_a0i_doFill(0, 0, true);
      std::get<0>(u0i_a0i_doFill) = LU.second.VectorIndex(0, 0, i);
      if (!matrix.IsZero(0, i))
      {
        std::get<1>(u0i_a0i_doFill) = matrix.VectorIndex(0, 0, i);
        std::get<2>(u0i_a0i_doFill) = false;
      }
      u0i_a0i_doFill_.push_back(u0i_a0i_doFill);
    }
    for (std::size_t k = 1; k < n; ++k)
    {
      std::tuple<std::size_t, std::size_t, std::size_t, char> nujk_nljk_u0k_doFill(0, 0, 0, true);
      bool has_u0k = !LU.second.IsZero(0, k);
      if (has_u0k)
      {
        std::get<2>(nujk_nljk_u0k_doFill) = LU.second.VectorIndex(0, 0, k);
        std::get<3>(nujk_nljk_u0k_doFill) = false;
      }
      for (std::size_t j = 1; j <= k; ++j)
      {
        if (LU.second.IsZero(j, k))
          continue;
        std::tuple<std::size_t, std::size_t, std::size_t, char, char> ujk_ajk_lj0_doFill_skipLU(0, 0, 0, true, true);
        std::get<0>(ujk_ajk_lj0_doFill_skipLU) = LU.second.VectorIndex(0, j, k);
        bool has_lj0_u0k = (!LU.first.IsZero(j, 0)) && has_u0k;
        if (!matrix.IsZero(j, k))
        {
          std::get<1>(ujk_ajk_lj0_doFill_skipLU) = matrix.VectorIndex(0, j, k);
          std::get<3>(ujk_ajk_lj0_doFill_skipLU) = false;
        }
        if (has_lj0_u0k)
        {
          std::get<2>(ujk_ajk_lj0_doFill_skipLU) = LU.first.VectorIndex(0, j, 0);
          std::get<4>(ujk_ajk_lj0_doFill_skipLU) = false;
        }
        ujk_ajk_lj0_doFill_skipLU_.push_back(ujk_ajk_lj0_doFill_skipLU);
        ++(std::get<0>(nujk_nljk_u0k_doFill));
      }
      for (std::size_t j = k + 1; j < n; ++j)
      {
        if (LU.first.IsZero(j, k))
          continue;
        std::tuple<std::size_t, std::size_t, std::size_t, char, char> ljk_ajk_lj0_doFill_skipLU(0, 0, 0, true, true);
        std::get<0>(ljk_ajk_lj0_doFill_skipLU) = LU.first.VectorIndex(0, j, k);
        bool has_lj0_u0k = (!LU.first.IsZero(j, 0)) && has_u0k;
        if (!matrix.IsZero(j, k))
        {
          std::get<1>(ljk_ajk_lj0_doFill_skipLU) = matrix.VectorIndex(0, j, k);
          std::get<3>(ljk_ajk_lj0_doFill_skipLU) = false;
        }
        if (has_lj0_u0k)
        {
          std::get<2>(ljk_ajk_lj0_doFill_skipLU) = LU.first.VectorIndex(0, j, 0);
          std::get<4>(ljk_ajk_lj0_doFill_skipLU) = false;
        }
        ljk_ajk_lj0_doFill_skipLU_.push_back(ljk_ajk_lj0_doFill_skipLU);
        ++(std::get<1>(nujk_nljk_u0k_doFill));
      }
      nujk_nljk_u0k_doFill_.push_back(nujk_nljk_u0k_doFill);
    }
    for (std::size_t i = 1; i < n; ++i)
    {
      std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> lii_uii_nj_nk(0, 0, 0, 0);
      std::get<0>(lii_uii_nj_nk) = LU.first.VectorIndex(0, i, i);
      std::get<1>(lii_uii_nj_nk) = LU.second.VectorIndex(0, i, i);
      // middle j loop to set L[j][i]
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (LU.first.IsZero(j, i))
          continue;
        lji_.push_back(LU.first.VectorIndex(0, j, i));
        ++(std::get<2>(lii_uii_nj_nk));
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
        ++(std::get<3>(lii_uii_nj_nk));
      }
      lii_uii_nj_nk_.push_back(lii_uii_nj_nk);
    }
  }

  template<class SparseMatrixPolicy>
  requires(
      SparseMatrixConcept<SparseMatrixPolicy>) inline std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LuDecompositionMozart::
      GetLUMatrices(const SparseMatrixPolicy& A, typename SparseMatrixPolicy::value_type initial_value)
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
          if(std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
            U_ids.insert(std::make_pair(j, k));
        for (std::size_t j = k + 1; j < n; ++j)
          if(std::find(L_ids.begin(), L_ids.end(), std::make_pair(j, i)) != L_ids.end())
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
  requires(!VectorizableSparse<SparseMatrixPolicy>) inline void LuDecompositionMozart::Decompose(
      const SparseMatrixPolicy& A,
      SparseMatrixPolicy& L,
      SparseMatrixPolicy& U) const
  {
    MICM_PROFILE_FUNCTION();
    const std::size_t n = A.NumRows();

    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto L_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto U_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto ujk_ajk_lj0_doFill_skipLU = ujk_ajk_lj0_doFill_skipLU_.begin();
      auto ljk_ajk_lj0_doFill_skipLU = ljk_ajk_lj0_doFill_skipLU_.begin();
      auto lii_uii_nj_nk = lii_uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();
      double Uii_inverse;
      double Uik;

      L_vector[std::get<0>(l00_u00_a00_)] = 1.0;
      U_vector[std::get<1>(l00_u00_a00_)] = A_vector[std::get<2>(l00_u00_a00_)];
      Uii_inverse = 1.0 / U_vector[std::get<1>(l00_u00_a00_)];      
      for (auto& li0_ai0_doFill : li0_ai0_doFill_)
      {
        if (std::get<2>(li0_ai0_doFill))
        {
          L_vector[std::get<0>(li0_ai0_doFill)] = 0;
        } else {
          L_vector[std::get<0>(li0_ai0_doFill)] = A_vector[std::get<1>(li0_ai0_doFill)] * Uii_inverse;
        }
      }
      for (auto& u0i_a0i_doFill : u0i_a0i_doFill_)
      {
        if (std::get<2>(u0i_a0i_doFill))
        {
          U_vector[std::get<0>(u0i_a0i_doFill)] = 0;
        } else {
          U_vector[std::get<0>(u0i_a0i_doFill)] = A_vector[std::get<1>(u0i_a0i_doFill)];
        }
      }
      for (auto& nujk_nljk_u0k_doFill : nujk_nljk_u0k_doFill_)
      {
        if (std::get<3>(nujk_nljk_u0k_doFill))
        {
          Uik = 0;
        }
        else
        {
          Uik = U_vector[std::get<2>(nujk_nljk_u0k_doFill)];
        }
        for (std::size_t j = 0; j < std::get<0>(nujk_nljk_u0k_doFill); ++j)
        {
          const std::size_t ujk = std::get<0>(*ujk_ajk_lj0_doFill_skipLU);
          if (std::get<3>(*ujk_ajk_lj0_doFill_skipLU))
          {
            U_vector[ujk] = 0;
          }
          else
          {
            U_vector[ujk] = A_vector[std::get<1>(*ujk_ajk_lj0_doFill_skipLU)];
          }
          if (!std::get<4>(*ujk_ajk_lj0_doFill_skipLU))
          {
            U_vector[ujk] -= L_vector[std::get<2>(*ujk_ajk_lj0_doFill_skipLU)] * Uik;
          }
          ++ujk_ajk_lj0_doFill_skipLU;
        }
        for (std::size_t j = 0; j < std::get<1>(nujk_nljk_u0k_doFill); ++j)
        {
          const std::size_t ljk = std::get<0>(*ljk_ajk_lj0_doFill_skipLU);
          if (std::get<3>(*ljk_ajk_lj0_doFill_skipLU))
          {
            L_vector[ljk] = 0;
          }
          else
          {
            L_vector[ljk] = A_vector[std::get<1>(*ljk_ajk_lj0_doFill_skipLU)];
          }
          if (!std::get<4>(*ljk_ajk_lj0_doFill_skipLU))
          {
            L_vector[ljk] -= L_vector[std::get<2>(*ljk_ajk_lj0_doFill_skipLU)] * Uik;
          }
          ++ljk_ajk_lj0_doFill_skipLU;
        }
      }
      for (std::size_t i = 1; i < n; ++i)
      {
        L_vector[std::get<0>(*lii_uii_nj_nk)] = 1.0;
        Uii_inverse = 1.0 / U_vector[std::get<1>(*lii_uii_nj_nk)];
        for (std::size_t j = 0; j < std::get<2>(*lii_uii_nj_nk); ++j)
        {
          L_vector[*lji] = L_vector[*lji] * Uii_inverse;
          ++lji;
        }
        for (std::size_t k = 0; k < std::get<3>(*lii_uii_nj_nk); ++k)
        {
          const double Uik = U_vector[std::get<2>(*nujk_nljk_uik)];
          for (std::size_t j = 0; j < std::get<0>(*nujk_nljk_uik); ++j)
          {
            U_vector[ujk_lji->first] -= L_vector[ujk_lji->second] * Uik;
            ++ujk_lji;
          }
          for (std::size_t j = 0 ; j < std::get<1>(*nujk_nljk_uik); ++j)
          {
            L_vector[ljk_lji->first] -= L_vector[ljk_lji->second] * Uik;
            ++ljk_lji;
          }
          ++nujk_nljk_uik;
        }
        ++lii_uii_nj_nk;
      }
    }
  }

  template<class SparseMatrixPolicy>
  requires(VectorizableSparse<SparseMatrixPolicy>) inline void LuDecompositionMozart::Decompose(
      const SparseMatrixPolicy& A,
      SparseMatrixPolicy& L,
      SparseMatrixPolicy& U) const
  {
    MICM_PROFILE_FUNCTION();

    const std::size_t n = A.NumRows();
    const std::size_t A_BlockSize = A.NumberOfBlocks();
    constexpr std::size_t A_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t A_GroupSizeOfFlatBlockSize = A.GroupSize(A.FlatBlockSize());
    const std::size_t L_GroupSizeOfFlatBlockSize = L.GroupSize(L.FlatBlockSize());
    const std::size_t U_GroupSizeOfFlatBlockSize = U.GroupSize(U.FlatBlockSize());
    double Uii_inverse[A_GroupVectorSize];
    double Uik[A_GroupVectorSize];

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A_BlockSize); ++i_group)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_group * A_GroupSizeOfFlatBlockSize);
      auto L_vector = std::next(L.AsVector().begin(), i_group * L_GroupSizeOfFlatBlockSize);
      auto U_vector = std::next(U.AsVector().begin(), i_group * U_GroupSizeOfFlatBlockSize);
      auto ujk_ajk_lj0_doFill_skipLU = ujk_ajk_lj0_doFill_skipLU_.begin();
      auto ljk_ajk_lj0_doFill_skipLU = ljk_ajk_lj0_doFill_skipLU_.begin();
      auto lii_uii_nj_nk = lii_uii_nj_nk_.begin();
      auto lji = lji_.begin();
      auto nujk_nljk_uik = nujk_nljk_uik_.begin();
      auto ujk_lji = ujk_lji_.begin();
      auto ljk_lji = ljk_lji_.begin();
      const std::size_t n_cells = std::min(A_GroupVectorSize, A_BlockSize - i_group * A_GroupVectorSize);
      
      for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
        L_vector[std::get<0>(l00_u00_a00_)+i_cell] = 1.0;
      for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
      {
        U_vector[std::get<1>(l00_u00_a00_)+i_cell] = A_vector[std::get<2>(l00_u00_a00_)+i_cell];
        Uii_inverse[i_cell] = 1.0 / U_vector[std::get<1>(l00_u00_a00_)+i_cell];
      }
      for (auto& li0_ai0_doFill : li0_ai0_doFill_)
      {
        if (std::get<2>(li0_ai0_doFill))
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[std::get<0>(li0_ai0_doFill)+i_cell] = 0;
        }
        else
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[std::get<0>(li0_ai0_doFill)+i_cell] = A_vector[std::get<1>(li0_ai0_doFill)+i_cell] * Uii_inverse[i_cell];
        }
      }
      for (auto& u0i_a0i_doFill : u0i_a0i_doFill_)
      {
        if (std::get<2>(u0i_a0i_doFill))
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            U_vector[std::get<0>(u0i_a0i_doFill)+i_cell] = 0;
        }
        else
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            U_vector[std::get<0>(u0i_a0i_doFill)+i_cell] = A_vector[std::get<1>(u0i_a0i_doFill)+i_cell];
        }
      }
      for (auto& nujk_nljk_u0k_doFill : nujk_nljk_u0k_doFill_)
      {
        if (std::get<3>(nujk_nljk_u0k_doFill))
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            Uik[i_cell] = 0;
        }
        else
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            Uik[i_cell] = U_vector[std::get<2>(nujk_nljk_u0k_doFill)+i_cell];
        }
        for (std::size_t j = 0; j < std::get<0>(nujk_nljk_u0k_doFill); ++j)
        {
          const std::size_t ujk = std::get<0>(*ujk_ajk_lj0_doFill_skipLU);
          if (std::get<3>(*ujk_ajk_lj0_doFill_skipLU))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[ujk+i_cell] = 0;
          }
          else
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[ujk+i_cell] = A_vector[std::get<1>(*ujk_ajk_lj0_doFill_skipLU)+i_cell];
          }
          if (!std::get<4>(*ujk_ajk_lj0_doFill_skipLU))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[ujk+i_cell] -= L_vector[std::get<2>(*ujk_ajk_lj0_doFill_skipLU)+i_cell] * Uik[i_cell];
          }
          ++ujk_ajk_lj0_doFill_skipLU;
        }
        for (std::size_t j = 0; j < std::get<1>(nujk_nljk_u0k_doFill); ++j)
        {
          const std::size_t ljk = std::get<0>(*ljk_ajk_lj0_doFill_skipLU);
          if (std::get<3>(*ljk_ajk_lj0_doFill_skipLU))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[ljk+i_cell] = 0;
          }
          else
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[ljk+i_cell] = A_vector[std::get<1>(*ljk_ajk_lj0_doFill_skipLU)+i_cell];
          }
          if (!std::get<4>(*ljk_ajk_lj0_doFill_skipLU))
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[ljk+i_cell] -= L_vector[std::get<2>(*ljk_ajk_lj0_doFill_skipLU)+i_cell] * Uik[i_cell];
          }
          ++ljk_ajk_lj0_doFill_skipLU;
        }
      }
      for (std::size_t i = 1; i < n; ++i)
      {
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          L_vector[std::get<0>(*lii_uii_nj_nk)+i_cell] = 1.0;
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          Uii_inverse[i_cell] = 1.0 / U_vector[std::get<1>(*lii_uii_nj_nk)+i_cell];
        for (std::size_t j = 0; j < std::get<2>(*lii_uii_nj_nk); ++j)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            L_vector[*lji+i_cell] = L_vector[*lji+i_cell] * Uii_inverse[i_cell];
          ++lji;
        }
        for (std::size_t k = 0; k < std::get<3>(*lii_uii_nj_nk); ++k)
        {
          const std::size_t uik = std::get<2>(*nujk_nljk_uik);
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            Uik[i_cell] = U_vector[uik+i_cell];
          for (std::size_t j = 0; j < std::get<0>(*nujk_nljk_uik); ++j)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              U_vector[ujk_lji->first+i_cell] -= L_vector[ujk_lji->second+i_cell] * Uik[i_cell];
            ++ujk_lji;
          }
          for (std::size_t j = 0; j < std::get<1>(*nujk_nljk_uik); ++j)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              L_vector[ljk_lji->first+i_cell] -= L_vector[ljk_lji->second+i_cell] * Uik[i_cell];
            ++ljk_lji;
          }
          ++nujk_nljk_uik;
        }
        ++lii_uii_nj_nk;
      }
    }
  }

}  // namespace micm
