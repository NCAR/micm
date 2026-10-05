// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart<SparseMatrixPolicy>::LuDecompositionMozart() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart<SparseMatrixPolicy>::LuDecompositionMozart(const SparseMatrixPolicy& matrix)
  {
    Initialize(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart<SparseMatrixPolicy>::LuDecompositionMozart(LuDecompositionMozart&& other) noexcept
      : lii_nuji_nlji_(std::move(other.lii_nuji_nlji_)),
        uji_aji_(std::move(other.uji_aji_)),
        lji_aji_(std::move(other.lji_aji_)),
        fill_uji_(std::move(other.fill_uji_)),
        fill_lji_(std::move(other.fill_lji_)),
        uii_nj_nk_(std::move(other.uii_nj_nk_)),
        lji_(std::move(other.lji_)),
        nujk_nljk_uik_(std::move(other.nujk_nljk_uik_)),
        ujk_lji_(std::move(other.ujk_lji_)),
        ljk_lji_(std::move(other.ljk_lji_)),
        views_(
            lii_nuji_nlji_,
            uji_aji_,
            lji_aji_,
            fill_uji_,
            fill_lji_,
            uii_nj_nk_,
            lji_,
            nujk_nljk_uik_,
            ujk_lji_,
            ljk_lji_)
  {
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart<SparseMatrixPolicy>& LuDecompositionMozart<SparseMatrixPolicy>::operator=(
      LuDecompositionMozart&& other) noexcept
  {
    if (this != &other)
    {
      lii_nuji_nlji_ = std::move(other.lii_nuji_nlji_);
      uji_aji_ = std::move(other.uji_aji_);
      lji_aji_ = std::move(other.lji_aji_);
      fill_uji_ = std::move(other.fill_uji_);
      fill_lji_ = std::move(other.fill_lji_);
      uii_nj_nk_ = std::move(other.uii_nj_nk_);
      lji_ = std::move(other.lji_);
      nujk_nljk_uik_ = std::move(other.nujk_nljk_uik_);
      ujk_lji_ = std::move(other.ujk_lji_);
      ljk_lji_ = std::move(other.ljk_lji_);
      views_ = Views(
          lii_nuji_nlji_, uji_aji_, lji_aji_, fill_uji_, fill_lji_, uii_nj_nk_, lji_, nujk_nljk_uik_, ujk_lji_, ljk_lji_);
    }
    return *this;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozart<SparseMatrixPolicy> LuDecompositionMozart<SparseMatrixPolicy>::Create(
      const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozart<SparseMatrixPolicy> lu_decomp{};
    lu_decomp.Initialize(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozart<SparseMatrixPolicy>::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    Index n = matrix.NumRows();
    auto LU = GetLUMatrices(matrix, initial_value, true);
    std::vector<IndexTrio> lii_nuji_nlji_temp;
    std::vector<IndexPair> uji_aji_temp;
    std::vector<IndexPair> lji_aji_temp;
    std::vector<Index> fill_uji_temp;
    std::vector<Index> fill_lji_temp;
    std::vector<IndexTrio> uii_nj_nk_temp;
    std::vector<Index> lji_temp;
    std::vector<IndexTrio> nujk_nljk_uik_temp;
    std::vector<IndexPair> ujk_lji_temp;
    std::vector<IndexPair> ljk_lji_temp;
    for (Index i = 0; i < n; ++i)
    {
      IndexTrio lii_nuji_nlji(0, 0, 0);
      lii_nuji_nlji.first_ = LU.first.VectorIndex(0, i, i);
      // set initial values for U matrix
      for (Index j = 0; j <= i; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!LU.second.IsZero(j, i))
          {
            fill_uji_temp.push_back(LU.second.VectorIndex(0, j, i));
          }
          continue;
        }
        uji_aji_temp.push_back({ LU.second.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i) });
        ++(lii_nuji_nlji.second_);
      }
      // set initial values for L matrix
      for (Index j = i + 1; j < n; ++j)
      {
        if (matrix.IsZero(j, i))
        {
          if (!LU.first.IsZero(j, i))
          {
            fill_lji_temp.push_back(LU.first.VectorIndex(0, j, i));
          }
          continue;
        }
        lji_aji_temp.push_back({ LU.first.VectorIndex(0, j, i), matrix.VectorIndex(0, j, i) });
        ++(lii_nuji_nlji.third_);
      }
      lii_nuji_nlji_temp.push_back(lii_nuji_nlji);
    }
    for (Index i = 0; i < matrix.NumRows(); ++i)
    {
      IndexTrio uii_nj_nk(0, 0, 0);
      uii_nj_nk.first_ = LU.second.VectorIndex(0, i, i);
      // middle j loop to set L[j][i]
      for (Index j = i + 1; j < n; ++j)
      {
        if (LU.first.IsZero(j, i))
        {
          continue;
        }
        lji_temp.push_back(LU.first.VectorIndex(0, j, i));
        ++(uii_nj_nk.second_);
      }
      // middle k loop to set U[j][k] and L[j][k]
      for (Index k = i + 1; k < n; ++k)
      {
        if (LU.second.IsZero(i, k))
        {
          continue;
        }
        IndexTrio nujk_nljk_uik(0, 0, 0);
        nujk_nljk_uik.third_ = LU.second.VectorIndex(0, i, k);
        // inner j loop to set U[j][k]
        for (Index j = i + 1; j <= k; ++j)
        {
          if (LU.first.IsZero(j, i))
          {
            continue;
          }
          IndexPair ujk_lji(0, 0);
          ujk_lji.first_ = LU.second.VectorIndex(0, j, k);
          ujk_lji.second_ = LU.first.VectorIndex(0, j, i);
          ujk_lji_temp.push_back(ujk_lji);
          ++(nujk_nljk_uik.first_);
        }
        // inner j loop to set L[j][k]
        for (Index j = k + 1; j < n; ++j)
        {
          if (LU.first.IsZero(j, i))
          {
            continue;
          }
          IndexPair ljk_lji(0, 0);
          ljk_lji.first_ = LU.first.VectorIndex(0, j, k);
          ljk_lji.second_ = LU.first.VectorIndex(0, j, i);
          ljk_lji_temp.push_back(ljk_lji);
          ++(nujk_nljk_uik.second_);
        }
        nujk_nljk_uik_temp.push_back(nujk_nljk_uik);
        ++(uii_nj_nk.third_);
      }
      uii_nj_nk_temp.push_back(uii_nj_nk);
    }
    lii_nuji_nlji_ = lii_nuji_nlji_temp;
    uji_aji_ = uji_aji_temp;
    lji_aji_ = lji_aji_temp;
    fill_uji_ = fill_uji_temp;
    fill_lji_ = fill_lji_temp;
    uii_nj_nk_ = uii_nj_nk_temp;
    lji_ = lji_temp;
    nujk_nljk_uik_ = nujk_nljk_uik_temp;
    ujk_lji_ = ujk_lji_temp;
    ljk_lji_ = ljk_lji_temp;
    lii_nuji_nlji_.CopyToDevice();
    uji_aji_.CopyToDevice();
    lji_aji_.CopyToDevice();
    fill_uji_.CopyToDevice();
    fill_lji_.CopyToDevice();
    uii_nj_nk_.CopyToDevice();
    lji_.CopyToDevice();
    nujk_nljk_uik_.CopyToDevice();
    ujk_lji_.CopyToDevice();
    ljk_lji_.CopyToDevice();
    views_ = Views(
        lii_nuji_nlji_, uji_aji_, lji_aji_, fill_uji_, fill_lji_, uii_nj_nk_, lji_, nujk_nljk_uik_, ujk_lji_, ljk_lji_);
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LuDecompositionMozart<SparseMatrixPolicy>::GetLUMatrices(
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
    auto L_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : L_ids)
    {
      L_builder = L_builder.WithElement(pair.first, pair.second);
    }
    auto U_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : U_ids)
    {
      U_builder = U_builder.WithElement(pair.first, pair.second);
    }
    std::pair<SparseMatrixPolicy, SparseMatrixPolicy> LU(
        SparseMatrixPolicy(L_builder, indexing_only), SparseMatrixPolicy(U_builder, indexing_only));
    return LU;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozart<SparseMatrixPolicy>::Decompose(
      const SparseMatrixPolicy& A,
      SparseMatrixPolicy& L,
      SparseMatrixPolicy& U) const
  {
    const Index n = A.NumRows();
    const auto& views = views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename SparseMatrix::ConstViewType& A_view,
            const typename SparseMatrix::ViewType& lower_view,
            const typename SparseMatrix::ViewType& upper_view) {
          auto uji_aji = views.uji_aji_.begin();
          auto lji_aji = views.lji_aji_.begin();
          auto uii_nj_nk = views.uii_nj_nk_.begin();
          auto lji = views.lji_.begin();
          auto nujk_nljk_uik = views.nujk_nljk_uik_.begin();
          auto ujk_lji = views.ujk_lji_.begin();
          auto ljk_lji = views.ljk_lji_.begin();
          auto Uii_inverse = A_view.GetBlockVariable();
          for (const auto& lii_nuji_nlji : views.lii_nuji_nlji_)
          {
            for (Index i = 0; i < lii_nuji_nlji.second_; ++i)
            {
              upper_view.Copy(upper_view.GetBlockView(uji_aji->first_), A_view.GetConstBlockView(uji_aji->second_));
              ++uji_aji;
            }
            lower_view.Fill(lower_view.GetBlockView(lii_nuji_nlji.first_), 1.0);
            for (Index i = 0; i < lii_nuji_nlji.third_; ++i)
            {
              lower_view.Copy(lower_view.GetBlockView(lji_aji->first_), A_view.GetConstBlockView(lji_aji->second_));
              ++lji_aji;
            }
          }
          for (const auto& fill_uji : views.fill_uji_)
          {
            upper_view.Fill(upper_view.GetBlockView(fill_uji), 0.0);
          }
          for (const auto& fill_lji : views.fill_lji_)
          {
            lower_view.Fill(lower_view.GetBlockView(fill_lji), 0.0);
          }
          for (Index i = 0; i < n; ++i)
          {
            upper_view.ForEachBlock(
                [](Real& uii_inv, const Real& uii) { uii_inv = 1.0 / uii; },
                Uii_inverse,
                upper_view.GetConstBlockView(uii_nj_nk->first_));
            for (Index ij = 0; ij < uii_nj_nk->second_; ++ij)
            {
              lower_view.ForEachBlock(
                  [](Real& lji, const Real& uii_inv) { lji *= uii_inv; }, lower_view.GetBlockView(*(lji++)), Uii_inverse);
            }
            for (Index ik = 0; ik < uii_nj_nk->third_; ++ik)
            {
              auto uik_view = upper_view.GetConstBlockView(nujk_nljk_uik->third_);
              for (Index ij = 0; ij < nujk_nljk_uik->first_; ++ij)
              {
                upper_view.ForEachBlock(
                    [](Real& ujk, const Real& lji, const Real& uik) { ujk -= lji * uik; },
                    upper_view.GetBlockView(ujk_lji->first_),
                    lower_view.GetConstBlockView(ujk_lji->second_),
                    uik_view);
                ++ujk_lji;
              }
              for (Index ij = 0; ij < nujk_nljk_uik->second_; ++ij)
              {
                lower_view.ForEachBlock(
                    [](Real& ljk, const Real& lji, const Real& uik) { ljk -= lji * uik; },
                    lower_view.GetBlockView(ljk_lji->first_),
                    lower_view.GetConstBlockView(ljk_lji->second_),
                    uik_view);
                ++ljk_lji;
              }
              ++nujk_nljk_uik;
            }
            ++uii_nj_nk;
          }
        },
        A,
        L,
        U)(A, L, U);
  }
}  // namespace micm
