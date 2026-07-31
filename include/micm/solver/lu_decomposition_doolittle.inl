// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  inline LuDecompositionDoolittle::LuDecompositionDoolittle() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle::LuDecompositionDoolittle(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle LuDecompositionDoolittle::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionDoolittle lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle::FillPattern LuDecompositionDoolittle::ComputeFillPattern(const SparseMatrixPolicy& A)
  {
    Index n = A.NumRows();
    FillPattern fp;
    fp.Arow_.assign(n, {});
    fp.Acol_.assign(n, {});
    fp.Lrow_.assign(n, {});
    fp.Urow_.assign(n, {});
    fp.Lcol_.assign(n, {});

    // Non-zero structure of the (sparse) input matrix A. Its rows are short, so the
    // O(n^2) IsZero scan here costs O(n * nnz(A)) -- negligible next to the dense
    // O(n^3) loop this whole routine replaces.
    for (Index r = 0; r < n; ++r)
    {
      for (Index c = 0; c < n; ++c)
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
    std::vector<std::vector<Index>> Ucol(n);  // Ucol[i] = rows j<i where U[j][i] != 0, ascending
    std::vector<Bool> seen(n, 0);
    std::vector<Index> touched;
    auto mark = [&](Index k)
    {
      if (!seen[k])
      {
        seen[k] = 1;
        touched.push_back(k);
      }
    };
    for (Index i = 0; i < n; ++i)
    {
      // Upper triangular matrix: columns k >= i that are non-zero in U row i.
      touched.clear();
      mark(i);  // unit/diagonal entry U[i][i]
      for (Index k : fp.Arow_[i])
      {
        if (k >= i)
        {
          mark(k);
        }
      }
      for (Index j : fp.Lrow_[i])  // j < i, L[i][j] != 0
      {
        for (Index k : fp.Urow_[j])  // k >= j, U[j][k] != 0
        {
          if (k >= i)
          {
            mark(k);
          }
        }
      }
      std::sort(touched.begin(), touched.end());
      for (Index k : touched)
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
      for (Index k : fp.Acol_[i])
      {
        if (k > i)
        {
          mark(k);
        }
      }
      for (Index j : Ucol[i])  // j < i, U[j][i] != 0
      {
        for (Index k : fp.Lcol_[j])  // k > j, L[k][j] != 0
        {
          if (k > i)
          {
            mark(k);
          }
        }
      }
      std::sort(touched.begin(), touched.end());
      for (Index k : touched)
      {
        fp.L_ids_.insert(std::make_pair(k, i));
        fp.Lcol_[i].push_back(k);
        fp.Lrow_[k].push_back(i);
        seen[k] = 0;
      }
    }
    return fp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    Index n = matrix.NumRows();
    FillPattern fp = ComputeFillPattern(matrix);
    // Build the (indexing-only) L and U matrices from the fill pattern so we can map
    // (row, column) positions to data-vector indices below.
    auto L_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(matrix.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.L_ids_)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(matrix.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.U_ids_)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LU(
        SparseMatrixPolicy(L_builder, true), SparseMatrixPolicy(U_builder, true));

    // O(1)-amortized membership on the sorted adjacency rows; bounded by the
    // factorization's own operation count (times a log factor) rather than O(n^3).
    auto contains = [](const std::vector<Index>& v, Index x) { return std::binary_search(v.begin(), v.end(), x); };

    for (Index i = 0; i < n; ++i)
    {
      std::pair<Index, Index> iLU(0, 0);
      // Upper triangular matrix: iterate only the non-zero columns of U row i.
      for (Index k : fp.Urow_[i])
      {
        Index nkj = 0;
        // j < i with L[i][j] != 0 and U[j][k] != 0, in ascending j order.
        for (Index j : fp.Lrow_[i])
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
      for (Index k : fp.Lcol_[i])
      {
        Index nkj = 0;
        // j < i with L[k][j] != 0 and U[j][i] != 0, in ascending j order. Lrow_[k] is
        // sorted, so stop once j reaches i.
        for (Index j : fp.Lrow_[k])
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

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LuDecompositionDoolittle::GetLUMatrices(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    Index n = A.NumRows();
    FillPattern fp = ComputeFillPattern(A);
    auto L_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.L_ids_)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : fp.U_ids_)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LU(
        SparseMatrixPolicy(L_builder, indexing_only), SparseMatrixPolicy(U_builder, indexing_only));
    return LU;
  }

  template<class SparseMatrixPolicy>
  inline void LuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    SparseMatrixPolicy::Function(
        [this](const auto&& A_view, auto&& lower_view, auto&& upper_view)
        {
          auto do_aik = do_aik_.begin();
          auto aik = aik_.begin();
          auto uik_nkj = uik_nkj_.begin();
          auto lij_ujk = lij_ujk_.begin();
          auto do_aki = do_aki_.begin();
          auto aki = aki_.begin();
          auto lki_nkj = lki_nkj_.begin();
          auto lkj_uji = lkj_uji_.begin();
          auto uii = uii_.begin();
          for (const auto& niLU : niLU_)
          {
            // Upper triangular matrix
            for (Index iU = 0; iU < niLU.second; ++iU)
            {
              auto uik_view = upper_view.GetBlockView(uik_nkj->first);
              if (*(do_aik++))
              {
                upper_view.Copy(uik_view, A_view.GetConstBlockView(*(aik++)));
              }
              else
              {
                upper_view.Fill(uik_view, 0.0);
              }
              for (Index ikj = 0; ikj < uik_nkj->second; ++ikj)
              {
                upper_view.ForEachBlock(
                    [](Real& uik, const Real& lij, const Real& ujk) { uik -= lij * ujk; },
                    uik_view,
                    lower_view.GetConstBlockView(lij_ujk->first),
                    upper_view.GetConstBlockView(lij_ujk->second));
                ++lij_ujk;
              }
              ++uik_nkj;
            }
            // Lower triangular matrix
            lower_view.Fill(lower_view.GetBlockView((lki_nkj++)->first), 1.0);
            for (Index iL = 0; iL < niLU.first; ++iL)
            {
              auto lki_view = lower_view.GetBlockView(lki_nkj->first);
              if (*(do_aki++))
              {
                lower_view.Copy(lki_view, A_view.GetConstBlockView(*(aki++)));
              }
              else
              {
                lower_view.Fill(lki_view, 0.0);
              }
              for (Index ikj = 0; ikj < lki_nkj->second; ++ikj)
              {
                lower_view.ForEachBlock(
                    [](Real& lki, const Real& lkj, const Real& uji) { lki -= lkj * uji; },
                    lki_view,
                    lower_view.GetConstBlockView(lkj_uji->first),
                    upper_view.GetConstBlockView(lkj_uji->second));
                ++lkj_uji;
              }
              lower_view.ForEachBlock(
                  [](Real& lki, const Real& uii) { lki /= uii; }, lki_view, upper_view.GetConstBlockView(*(uii++)));
              ++lki_nkj;
            }
          }
        },
        A,
        L,
        U)(A, L, U);
  }
}  // namespace micm
