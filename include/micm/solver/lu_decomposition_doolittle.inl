// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionDoolittle::LuDecompositionDoolittle() = default;

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

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle::FillPattern LuDecompositionDoolittle::ComputeFillPattern(const SparseMatrixPolicy& A)
  {
    std::size_t n = A.NumRows();
    FillPattern fp;
    fp.Arow_.assign(n, {});
    fp.Acol_.assign(n, {});
    fp.Lrow_.assign(n, {});
    fp.Urow_.assign(n, {});
    fp.Lcol_.assign(n, {});

    // Non-zero structure of the (sparse) input matrix A. Its rows are short, so the
    // O(n^2) IsZero scan here costs O(n * nnz(A)) -- negligible next to the dense
    // O(n^3) loop this whole routine replaces.
    for (std::size_t r = 0; r < n; ++r)
    {
      for (std::size_t c = 0; c < n; ++c)
      {
        if (!A.IsZero(r, c))
        {
          fp.Arow_[r].push_back(c);
          fp.Acol_[c].push_back(r);
        }
      }
    }

    // Symbolic factorization processed one row at a time in increasing-i order. By
    // the time row i is reached, L[i][.] (j<i), U[j][.] and L[.][j] (j<i) are known,
    // so the fill of U row i and L column i is found by walking only the relevant
    // non-zeros: U[i][k] fills iff A[i][k] != 0, k == i, or some j<i has L[i][j] and
    // U[j][k]; L[k][i] fills iff A[k][i] != 0 or some j<i has L[k][j] and U[j][i].
    std::vector<std::vector<std::size_t>> Ucol(n);  // Ucol[i] = rows j<i where U[j][i] != 0, ascending
    std::vector<char> seen(n, 0);
    std::vector<std::size_t> touched;
    auto mark = [&](std::size_t k)
    {
      if (!seen[k])
      {
        seen[k] = 1;
        touched.push_back(k);
      }
    };
    for (std::size_t i = 0; i < n; ++i)
    {
      // Upper triangular matrix: columns k >= i that are non-zero in U row i.
      touched.clear();
      mark(i);  // unit/diagonal entry U[i][i]
      for (std::size_t k : fp.Arow_[i])
      {
        if (k >= i)
        {
          mark(k);
        }
      }
      for (std::size_t j : fp.Lrow_[i])  // j < i, L[i][j] != 0
      {
        for (std::size_t k : fp.Urow_[j])  // k >= j, U[j][k] != 0
        {
          if (k >= i)
          {
            mark(k);
          }
        }
      }
      std::sort(touched.begin(), touched.end());
      for (std::size_t k : touched)
      {
        fp.U_ids_.insert(std::make_pair(i, k));
        fp.Urow_[i].push_back(k);
        Ucol[k].push_back(i);
        seen[k] = 0;
      }

      // Lower triangular matrix: rows k > i non-zero in L column i, plus the unit
      // diagonal L[i][i].
      fp.L_ids_.insert(std::make_pair(i, i));
      touched.clear();
      for (std::size_t k : fp.Acol_[i])
      {
        if (k > i)
        {
          mark(k);
        }
      }
      for (std::size_t j : Ucol[i])  // j < i, U[j][i] != 0
      {
        for (std::size_t k : fp.Lcol_[j])  // k > j, L[k][j] != 0
        {
          if (k > i)
          {
            mark(k);
          }
        }
      }
      std::sort(touched.begin(), touched.end());
      for (std::size_t k : touched)
      {
        fp.L_ids_.insert(std::make_pair(k, i));
        fp.Lcol_[i].push_back(k);
        fp.Lrow_[k].push_back(i);
        seen[k] = 0;
      }
    }
    return fp;
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    std::size_t n = matrix.NumRows();
    FillPattern fp = ComputeFillPattern(matrix);
    // Build the (indexing-only) L and U matrices from the fill pattern so we can map
    // (row, column) positions to data-vector indices below.
    auto L_builder = LMatrixPolicy::Create(n).SetNumberOfBlocks(matrix.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.L_ids_)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = UMatrixPolicy::Create(n).SetNumberOfBlocks(matrix.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.U_ids_)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<LMatrixPolicy, UMatrixPolicy> LU(LMatrixPolicy(L_builder, true), UMatrixPolicy(U_builder, true));

    // O(1)-amortized membership on the sorted adjacency rows; bounded by the
    // factorization's own operation count (times a log factor) rather than O(n^3).
    auto contains = [](const std::vector<std::size_t>& v, std::size_t x)
    { return std::binary_search(v.begin(), v.end(), x); };

    for (std::size_t i = 0; i < n; ++i)
    {
      std::pair<std::size_t, std::size_t> iLU(0, 0);
      // Upper triangular matrix: iterate only the non-zero columns of U row i.
      for (std::size_t k : fp.Urow_[i])
      {
        std::size_t nkj = 0;
        // j < i with L[i][j] != 0 and U[j][k] != 0, in ascending j order.
        for (std::size_t j : fp.Lrow_[i])
        {
          if (!contains(fp.Urow_[j], k))
          {
            continue;
          }
          ++nkj;
          lij_ujk_.push_back(std::make_pair(LU.first.VectorIndex(0, i, j), LU.second.VectorIndex(0, j, k)));
        }
        if (contains(fp.Arow_[i], k))
        {
          do_aik_.push_back(true);
          aik_.push_back(matrix.VectorIndex(0, i, k));
        }
        else
        {
          do_aik_.push_back(false);
        }
        uik_nkj_.push_back(std::make_pair(LU.second.VectorIndex(0, i, k), nkj));
        ++(iLU.second);
      }
      // Lower triangular matrix: iterate only the non-zero rows of L column i.
      lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, i, i), 0));
      for (std::size_t k : fp.Lcol_[i])
      {
        std::size_t nkj = 0;
        // j < i with L[k][j] != 0 and U[j][i] != 0, in ascending j order. Lrow_[k] is
        // sorted, so stop once j reaches i.
        for (std::size_t j : fp.Lrow_[k])
        {
          if (j >= i)
          {
            break;
          }
          if (!contains(fp.Urow_[j], i))
          {
            continue;
          }
          ++nkj;
          lkj_uji_.push_back(std::make_pair(LU.first.VectorIndex(0, k, j), LU.second.VectorIndex(0, j, i)));
        }
        if (contains(fp.Acol_[i], k))
        {
          do_aki_.push_back(true);
          aki_.push_back(matrix.VectorIndex(0, k, i));
        }
        else
        {
          do_aki_.push_back(false);
        }
        uii_.push_back(LU.second.VectorIndex(0, i, i));
        lki_nkj_.push_back(std::make_pair(LU.first.VectorIndex(0, k, i), nkj));
        ++(iLU.first);
      }
      niLU_.push_back(iLU);
    }
    uii_.push_back(LU.second.VectorIndex(0, n - 1, n - 1));
  }

  template<class SparseMatrixPolicy, class LMatrixPolicy, class UMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<LMatrixPolicy, UMatrixPolicy> LuDecompositionDoolittle::GetLUMatrices(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    std::size_t n = A.NumRows();
    FillPattern fp = ComputeFillPattern(A);
    auto L_builder = LMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.L_ids_)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = UMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.U_ids_)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<LMatrixPolicy, UMatrixPolicy> LU(
        LMatrixPolicy(L_builder, indexing_only), UMatrixPolicy(U_builder, indexing_only));
    return LU;
  }

  template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    // Loop over blocks
    for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
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
      for (const auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          if (*(do_aik++))
          {
            U_vector[uik_nkj->first] = A_vector[*(aik++)];
          }
          else
          {
            U_vector[uik_nkj->first] = 0;
          }
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
          {
            L_vector[lki_nkj->first] = A_vector[*(aki++)];
          }
          else
          {
            L_vector[lki_nkj->first] = 0;
          }
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

  template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    const std::size_t A_BlockSize = A.NumberOfBlocks();
    constexpr std::size_t A_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t A_GroupSizeOfFlatBlockSize = A.GroupSize();
    const std::size_t L_GroupSizeOfFlatBlockSize = L.GroupSize();
    const std::size_t U_GroupSizeOfFlatBlockSize = U.GroupSize();

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < A.NumberOfGroups(A_BlockSize); ++i_group)
    {
      auto A_vector = std::next(A.AsVector().begin(), i_group * A_GroupSizeOfFlatBlockSize);
      auto L_vector = std::next(L.AsVector().begin(), i_group * L_GroupSizeOfFlatBlockSize);
      auto U_vector = std::next(U.AsVector().begin(), i_group * U_GroupSizeOfFlatBlockSize);
      auto do_aik = do_aik_.begin();
      auto aik = aik_.begin();
      auto uik_nkj = uik_nkj_.begin();
      auto lij_ujk = lij_ujk_.begin();
      auto do_aki = do_aki_.begin();
      auto aki = aki_.begin();
      auto lki_nkj = lki_nkj_.begin();
      auto lkj_uji = lkj_uji_.begin();
      auto uii = uii_.begin();
      const std::size_t n_cells = std::min(A_GroupVectorSize, A_BlockSize - i_group * A_GroupVectorSize);
      for (const auto& inLU : niLU_)
      {
        // Upper trianglur matrix
        for (std::size_t iU = 0; iU < inLU.second; ++iU)
        {
          const std::size_t uik_nkj_first = uik_nkj->first;
          if (*(do_aik++))
          {
            std::copy(A_vector + *aik, A_vector + *aik + n_cells, U_vector + uik_nkj_first);
            ++aik;
          }
          else
          {
            std::fill(U_vector + uik_nkj_first, U_vector + uik_nkj_first + n_cells, 0);
          }
          for (std::size_t ikj = 0; ikj < uik_nkj->second; ++ikj)
          {
            const std::size_t lij_ujk_first = lij_ujk->first;
            const std::size_t lij_ujk_second = lij_ujk->second;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              U_vector[uik_nkj_first + i_cell] -= L_vector[lij_ujk_first + i_cell] * U_vector[lij_ujk_second + i_cell];
            }
            ++lij_ujk;
          }
          ++uik_nkj;
        }
        // Lower triangular matrix
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          L_vector[lki_nkj->first + i_cell] = 1.0;
        }
        ++lki_nkj;
        for (std::size_t iL = 0; iL < inLU.first; ++iL)
        {
          const std::size_t lki_nkj_first = lki_nkj->first;
          if (*(do_aki++))
          {
            std::copy(A_vector + *aki, A_vector + *aki + n_cells, L_vector + lki_nkj_first);
            ++aki;
          }
          else
          {
            std::fill(L_vector + lki_nkj_first, L_vector + lki_nkj_first + n_cells, 0);
          }
          for (std::size_t ikj = 0; ikj < lki_nkj->second; ++ikj)
          {
            const std::size_t lkj_uji_first = lkj_uji->first;
            const std::size_t lkj_uji_second = lkj_uji->second;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              L_vector[lki_nkj_first + i_cell] -= L_vector[lkj_uji_first + i_cell] * U_vector[lkj_uji_second + i_cell];
            }
            ++lkj_uji;
          }
          const std::size_t uii_deref = *uii;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            L_vector[lki_nkj_first + i_cell] /= U_vector[uii_deref + i_cell];
          }
          ++lki_nkj;
          ++uii;
        }
      }
    }
  }

}  // namespace micm
