// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionDoolittle::LuDecompositionDoolittle()
  {
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle::LuDecompositionDoolittle(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle LuDecompositionDoolittle::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionDoolittle lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(
        matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    std::size_t n = matrix.NumRows();
    auto lu = GetLUMatrices<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, initial_value, true);
    for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    {
      std::pair<std::size_t, std::size_t> iLU(0, 0);
      // Upper triangular matrix
      for (std::size_t k = i; k < n; ++k)
      {
        std::size_t nkj = 0;
        for (std::size_t j = 0; j < i; ++j)
        {
          if (lu.first.IsZero(i, j) || lu.second.IsZero(j, k))
          {
            continue;
          }
          ++nkj;
          lij_ujk_.push_back(std::make_pair(lu.first.VectorIndex(0, i, j), lu.second.VectorIndex(0, j, k)));
        }
        if (matrix.IsZero(i, k))
        {
          if (nkj == 0 && k != i)
          {
            continue;
          }
          do_aik_.push_back(false);
        }
        else
        {
          do_aik_.push_back(true);
          aik_.push_back(matrix.VectorIndex(0, i, k));
        }
        uik_nkj_.push_back(std::make_pair(lu.second.VectorIndex(0, i, k), nkj));
        ++(iLU.second);
      }
      // Lower triangular matrix
      lki_nkj_.push_back(std::make_pair(lu.first.VectorIndex(0, i, i), 0));
      for (std::size_t k = i + 1; k < n; ++k)
      {
        std::size_t nkj = 0;
        for (std::size_t j = 0; j < i; ++j)
        {
          if (lu.first.IsZero(k, j) || lu.second.IsZero(j, i))
          {
            continue;
          }
          ++nkj;
          lkj_uji_.push_back(std::make_pair(lu.first.VectorIndex(0, k, j), lu.second.VectorIndex(0, j, i)));
        }
        if (matrix.IsZero(k, i))
        {
          if (nkj == 0)
          {
            continue;
          }
          do_aki_.push_back(false);
        }
        else
        {
          do_aki_.push_back(true);
          aki_.push_back(matrix.VectorIndex(0, k, i));
        }
        uii_.push_back(lu.second.VectorIndex(0, i, i));
        lki_nkj_.push_back(std::make_pair(lu.first.VectorIndex(0, k, i), nkj));
        ++(iLU.first);
      }
      niLU_.push_back(iLU);
    }
    uii_.push_back(lu.second.VectorIndex(0, n - 1, n - 1));
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<LMatrixPolicy, UMatrixPolicy> LuDecompositionDoolittle::GetLUMatrices(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    std::size_t n = A.NumRows();
    std::set<std::pair<std::size_t, std::size_t>> L_ids, U_ids;
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
  inline void LuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
    {
      auto a_vector = std::next(A.AsVector().begin(), i_block * A.FlatBlockSize());
      auto l_vector = std::next(L.AsVector().begin(), i_block * L.FlatBlockSize());
      auto u_vector = std::next(U.AsVector().begin(), i_block * U.FlatBlockSize());
      auto do_aik = do_aik_.begin();
      auto aik = aik_.begin();
      auto uik_nkj = uik_nkj_.begin();
      auto lij_ujk = lij_ujk_.begin();
      auto do_aki = do_aki_.begin();
      auto aki = aki_.begin();
      auto lki_nkj = lki_nkj_.begin();
      auto lkj_uji = lkj_uji_.begin();
      auto uii = uii_.begin();
      for (const auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          if (*(do_aik++))
          {
            u_vector[uik_nkj->first] = a_vector[*(aik++)];
          }
          else
          {
            u_vector[uik_nkj->first] = 0;
          }
          for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
          {
            u_vector[uik_nkj->first] -= l_vector[lij_ujk->first] * u_vector[lij_ujk->second];
            ++lij_ujk;
          }
          ++uik_nkj;
        }
        // Lower triangular matrix
        l_vector[(lki_nkj++)->first] = 1.0;
        for (std::size_t iL = 0; iL < inLU.first; ++iL)
        {
          if (*(do_aki++))
          {
            l_vector[lki_nkj->first] = a_vector[*(aki++)];
          }
          else
          {
            l_vector[lki_nkj->first] = 0;
          }
          for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
          {
            l_vector[lki_nkj->first] -= l_vector[lkj_uji->first] * u_vector[lkj_uji->second];
            ++lkj_uji;
          }
          l_vector[lki_nkj->first] /= u_vector[*uii];
          ++lki_nkj;
          ++uii;
        }
      }
    }
  }

  template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    const std::size_t A_BLOCK_SIZE = A.NumberOfBlocks();
    constexpr std::size_t A_GROUP_VECTOR_SIZE = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t A_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = A.GroupSize();
    const std::size_t L_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = L.GroupSize();
    const std::size_t U_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = U.GroupSize();

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A_BLOCK_SIZE); ++i_group)
    {
      auto a_vector = std::next(A.AsVector().begin(), i_group * A_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto l_vector = std::next(L.AsVector().begin(), i_group * L_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto u_vector = std::next(U.AsVector().begin(), i_group * U_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      auto do_aik = do_aik_.begin();
      auto aik = aik_.begin();
      auto uik_nkj = uik_nkj_.begin();
      auto lij_ujk = lij_ujk_.begin();
      auto do_aki = do_aki_.begin();
      auto aki = aki_.begin();
      auto lki_nkj = lki_nkj_.begin();
      auto lkj_uji = lkj_uji_.begin();
      auto uii = uii_.begin();
      const std::size_t N_CELLS = std::min(A_GROUP_VECTOR_SIZE, A_BLOCK_SIZE - i_group * A_GROUP_VECTOR_SIZE);
      for (const auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          const std::size_t uik_nkj_first = uik_nkj->first;
          if (*(do_aik++))
          {
            std::copy(a_vector + *aik, a_vector + *aik + N_CELLS, u_vector + uik_nkj_first);
            ++aik;
          }
          else
          {
            std::fill(u_vector + uik_nkj_first, u_vector + uik_nkj_first + N_CELLS, 0);
          }
          for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
          {
            const std::size_t lij_ujk_first = lij_ujk->first;
            const std::size_t lij_ujk_second = lij_ujk->second;
            for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
              u_vector[uik_nkj_first + i_cell] -= l_vector[lij_ujk_first + i_cell] * u_vector[lij_ujk_second + i_cell];
            ++lij_ujk;
          }
          ++uik_nkj;
        }
        // Lower triangular matrix
        for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          l_vector[lki_nkj->first + i_cell] = 1.0;
        ++lki_nkj;
        for (std::size_t iL = 0; iL < inLU.first; ++iL)
        {
          const std::size_t lki_nkj_first = lki_nkj->first;
          if (*(do_aki++))
          {
            std::copy(a_vector + *aki, a_vector + *aki + N_CELLS, l_vector + lki_nkj_first);
            ++aki;
          }
          else
          {
            std::fill(l_vector + lki_nkj_first, l_vector + lki_nkj_first + N_CELLS, 0);
          }
          for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
          {
            const std::size_t lkj_uji_first = lkj_uji->first;
            const std::size_t lkj_uji_second = lkj_uji->second;
            for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
              l_vector[lki_nkj_first + i_cell] -= l_vector[lkj_uji_first + i_cell] * u_vector[lkj_uji_second + i_cell];
            ++lkj_uji;
          }
          const std::size_t uii_deref = *uii;
          for (std::size_t i_cell = 0; i_cell < N_CELLS; ++i_cell)
          {
            l_vector[lki_nkj_first + i_cell] /= u_vector[uii_deref + i_cell];
          }
          ++lki_nkj;
          ++uii;
        }
      }
    }
  }

}  // namespace micm
