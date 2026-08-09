// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy>::LuDecompositionDoolittle() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy>::LuDecompositionDoolittle(const SparseMatrixPolicy& matrix)
  {
    Initialize(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy>::LuDecompositionDoolittle(LuDecompositionDoolittle&& other) noexcept
      : niLU_(std::move(other.niLU_)),
        do_aik_(std::move(other.do_aik_)),
        aik_(std::move(other.aik_)),
        uik_nkj_(std::move(other.uik_nkj_)),
        lij_ujk_(std::move(other.lij_ujk_)),
        do_aki_(std::move(other.do_aki_)),
        aki_(std::move(other.aki_)),
        lki_nkj_(std::move(other.lki_nkj_)),
        lkj_uji_(std::move(other.lkj_uji_)),
        uii_(std::move(other.uii_)),
        views_(niLU_, do_aik_, aik_, uik_nkj_, lij_ujk_, do_aki_, aki_, lki_nkj_, lkj_uji_, uii_)
  {
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy>& LuDecompositionDoolittle<SparseMatrixPolicy>::operator=(LuDecompositionDoolittle&& other) noexcept
  {
      if (this != &other)
      {
          niLU_ = std::move(other.niLU_);
          do_aik_ = std::move(other.do_aik_);
          aik_ = std::move(other.aik_);
          uik_nkj_ = std::move(other.uik_nkj_);
          lij_ujk_ = std::move(other.lij_ujk_);
          do_aki_ = std::move(other.do_aki_);
          aki_ = std::move(other.aki_);
          lki_nkj_ = std::move(other.lki_nkj_);
          lkj_uji_ = std::move(other.lkj_uji_);
          uii_ = std::move(other.uii_);
          views_       = Views(niLU_, do_aik_, aik_, uik_nkj_, lij_ujk_, do_aki_, aki_, lki_nkj_, lkj_uji_, uii_);
      }
      return *this;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy> LuDecompositionDoolittle<SparseMatrixPolicy>::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionDoolittle<SparseMatrixPolicy> lu_decomp{};
    lu_decomp.Initialize(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittle<SparseMatrixPolicy>::FillPattern LuDecompositionDoolittle<SparseMatrixPolicy>::ComputeFillPattern(const SparseMatrixPolicy& A)
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
  inline void LuDecompositionDoolittle<SparseMatrixPolicy>::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
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

    std::vector<IndexPair> niLU_temp;
    std::vector<Bool> do_aik_temp;
    std::vector<Index> aik_temp;
    std::vector<IndexPair> uik_nkj_temp;
    std::vector<IndexPair> lij_ujk_temp;
    std::vector<Bool> do_aki_temp;
    std::vector<Index> aki_temp;
    std::vector<IndexPair> lki_nkj_temp;
    std::vector<IndexPair> lkj_uji_temp;
    std::vector<Index> uii_temp; 

    for (Index i = 0; i < n; ++i)
    {
      IndexPair iLU(0, 0);
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
          lij_ujk_temp.push_back({LU.first.VectorIndex(0, i, j), LU.second.VectorIndex(0, j, k)});
        }
        if (contains(fp.Arow_[i], k))
        {
          do_aik_temp.push_back(true);
          aik_temp.push_back(matrix.VectorIndex(0, i, k));
        }
        else
        {
          do_aik_temp.push_back(false);
        }
        uik_nkj_temp.push_back({LU.second.VectorIndex(0, i, k), nkj});
        ++(iLU.second_);
      }
      // Lower triangular matrix: iterate only the non-zero rows of L column i.
      lki_nkj_temp.push_back({LU.first.VectorIndex(0, i, i), 0});
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
          lkj_uji_temp.push_back({LU.first.VectorIndex(0, k, j), LU.second.VectorIndex(0, j, i)});
        }
        if (contains(fp.Acol_[i], k))
        {
          do_aki_temp.push_back(true);
          aki_temp.push_back(matrix.VectorIndex(0, k, i));
        }
        else
        {
          do_aki_temp.push_back(false);
        }
        uii_temp.push_back(LU.second.VectorIndex(0, i, i));
        lki_nkj_temp.push_back({LU.first.VectorIndex(0, k, i), nkj});
        ++(iLU.first_);
      }
      niLU_temp.push_back(iLU);
    }
    uii_temp.push_back(LU.second.VectorIndex(0, n - 1, n - 1));
    niLU_ = niLU_temp;
    do_aik_ = do_aik_temp;
    aik_ = aik_temp;
    uik_nkj_ = uik_nkj_temp;
    lij_ujk_ = lij_ujk_temp;
    do_aki_ = do_aki_temp;
    aki_ = aki_temp;
    lki_nkj_ = lki_nkj_temp;
    lkj_uji_ = lkj_uji_temp;
    uii_ = uii_temp;
    niLU_.CopyToDevice();
    do_aik_.CopyToDevice();
    aik_.CopyToDevice();
    uik_nkj_.CopyToDevice();
    lij_ujk_.CopyToDevice();
    do_aki_.CopyToDevice();
    aki_.CopyToDevice();
    lki_nkj_.CopyToDevice();
    lkj_uji_.CopyToDevice();
    uii_.CopyToDevice();
    views_ = Views(niLU_, do_aik_, aik_, uik_nkj_, lij_ujk_, do_aki_, aki_, lki_nkj_, lkj_uji_, uii_);
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LuDecompositionDoolittle<SparseMatrixPolicy>::GetLUMatrices(
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
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittle<SparseMatrixPolicy>::Decompose(const SparseMatrixPolicy& A, SparseMatrixPolicy& L, SparseMatrixPolicy& U) const
  {
    const auto& views = views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(typename SparseMatrix::ConstViewType A_view, typename SparseMatrix::ViewType lower_view, typename SparseMatrix::ViewType upper_view)
        {
          auto do_aik = views.do_aik_.begin();
          auto aik = views.aik_.begin();
          auto uik_nkj = views.uik_nkj_.begin();
          auto lij_ujk = views.lij_ujk_.begin();
          auto do_aki = views.do_aki_.begin();
          auto aki = views.aki_.begin();
          auto lki_nkj = views.lki_nkj_.begin();
          auto lkj_uji = views.lkj_uji_.begin();
          auto uii = views.uii_.begin();
          for (const auto& niLU : views.niLU_)
          {
            // Upper triangular matrix
            for (Index iU = 0; iU < niLU.second_; ++iU)
            {
              auto uik_view = upper_view.GetBlockView(uik_nkj->first_);
              if (*(do_aik++))
              {
                upper_view.Copy(uik_view, A_view.GetConstBlockView(*(aik++)));
              }
              else
              {
                upper_view.Fill(uik_view, 0.0);
              }
              for (Index ikj = 0; ikj < uik_nkj->second_; ++ikj)
              {
                upper_view.ForEachBlock(
                    [](Real& uik, const Real& lij, const Real& ujk) { uik -= lij * ujk; },
                    uik_view,
                    lower_view.GetConstBlockView(lij_ujk->first_),
                    upper_view.GetConstBlockView(lij_ujk->second_));
                ++lij_ujk;
              }
              ++uik_nkj;
            }
            // Lower triangular matrix
            lower_view.Fill(lower_view.GetBlockView((lki_nkj++)->first_), 1.0);
            for (Index iL = 0; iL < niLU.first_; ++iL)
            {
              auto lki_view = lower_view.GetBlockView(lki_nkj->first_);
              if (*(do_aki++))
              {
                lower_view.Copy(lki_view, A_view.GetConstBlockView(*(aki++)));
              }
              else
              {
                lower_view.Fill(lki_view, 0.0);
              }
              for (Index ikj = 0; ikj < lki_nkj->second_; ++ikj)
              {
                lower_view.ForEachBlock(
                    [](Real& lki, const Real& lkj, const Real& uji) { lki -= lkj * uji; },
                    lki_view,
                    lower_view.GetConstBlockView(lkj_uji->first_),
                    upper_view.GetConstBlockView(lkj_uji->second_));
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
