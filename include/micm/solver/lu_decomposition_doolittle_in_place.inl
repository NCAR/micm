// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::LuDecompositionDoolittleInPlace() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::LuDecompositionDoolittleInPlace(
      const SparseMatrixPolicy& matrix)
  {
    Initialize(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::LuDecompositionDoolittleInPlace(
      LuDecompositionDoolittleInPlace&& other) noexcept
      : nik_nki_aii_(std::move(other.nik_nki_aii_)),
        aik_njk_(std::move(other.aik_njk_)),
        aij_ajk_(std::move(other.aij_ajk_)),
        aki_nji_(std::move(other.aki_nji_)),
        akj_aji_(std::move(other.akj_aji_)),
        views_(nik_nki_aii_, aik_njk_, aij_ajk_, aki_nji_, akj_aji_)
  {
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace<SparseMatrixPolicy>& LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::operator=(
      LuDecompositionDoolittleInPlace&& other) noexcept
  {
    if (this != &other)
    {
      nik_nki_aii_ = std::move(other.nik_nki_aii_);
      aik_njk_ = std::move(other.aik_njk_);
      aij_ajk_ = std::move(other.aij_ajk_);
      aki_nji_ = std::move(other.aki_nji_);
      akj_aji_ = std::move(other.akj_aji_);
      views_ = Views(nik_nki_aii_, aik_njk_, aij_ajk_, aki_nji_, akj_aji_);
    }
    return *this;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace<SparseMatrixPolicy> LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::Create(
      const SparseMatrixPolicy& matrix)
  {
    LuDecompositionDoolittleInPlace<SparseMatrixPolicy> lu_decomp{};
    lu_decomp.Initialize(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::Initialize(
      const SparseMatrixPolicy& matrix,
      auto initial_value)
  {
    Index n = matrix.NumRows();
    auto ALU = GetLUMatrix(matrix, initial_value, true);
    std::vector<IndexTrio> nik_nki_aii_temp;
    std::vector<IndexPair> aik_njk_temp;
    std::vector<IndexPair> aij_ajk_temp;
    std::vector<IndexPair> aki_nji_temp;
    std::vector<IndexPair> akj_aji_temp;
    for (Index i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      IndexTrio nik_nki_aii(0, 0, ALU.VectorIndex(0, i, i));
      for (Index k = i; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
        {
          continue;
        }
        IndexPair aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (Index j = 0; j < i; ++j)
        {
          if (ALU.IsZero(i, j) || ALU.IsZero(j, k))
          {
            continue;
          }
          aij_ajk_temp.push_back({ ALU.VectorIndex(0, i, j), ALU.VectorIndex(0, j, k) });
          ++(aik_njk.second_);
        }
        aik_njk_temp.push_back(aik_njk);
        ++(nik_nki_aii.first_);
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(k, i))
        {
          continue;
        }
        IndexPair aki_nji(ALU.VectorIndex(0, k, i), 0);
        for (Index j = 0; j < i; ++j)
        {
          if (ALU.IsZero(k, j) || ALU.IsZero(j, i))
          {
            continue;
          }
          akj_aji_temp.push_back({ ALU.VectorIndex(0, k, j), ALU.VectorIndex(0, j, i) });
          ++(aki_nji.second_);
        }
        aki_nji_temp.push_back(aki_nji);
        ++(nik_nki_aii.second_);
      }
      nik_nki_aii_temp.push_back(nik_nki_aii);
    }
    nik_nki_aii_ = nik_nki_aii_temp;
    aik_njk_ = aik_njk_temp;
    aij_ajk_ = aij_ajk_temp;
    aki_nji_ = aki_nji_temp;
    akj_aji_ = akj_aji_temp;
    nik_nki_aii_.CopyToDevice();
    aik_njk_.CopyToDevice();
    aij_ajk_.CopyToDevice();
    aki_nji_.CopyToDevice();
    akj_aji_.CopyToDevice();
    views_ = Views(nik_nki_aii_, aik_njk_, aij_ajk_, aki_nji_, akj_aji_);
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline SparseMatrixPolicy LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::GetLUMatrix(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    Index n = A.NumRows();
    std::set<std::pair<Index, Index>> ALU_ids;
    for (Index i = 0; i < n; ++i)
    {
      // Upper triangular matrix
      for (Index k = i; k < n; ++k)
      {
        if (!A.IsZero(i, k) || k == i)
        {
          ALU_ids.insert(std::make_pair(i, k));
          continue;
        }
        for (Index j = 0; j < i; ++j)
        {
          if (ALU_ids.contains(std::make_pair(i, j)) && ALU_ids.contains(std::make_pair(j, k)))
          {
            ALU_ids.insert(std::make_pair(i, k));
            break;
          }
        }
      }
      // Lower triangular matrix
      for (Index k = i; k < n; ++k)
      {
        if (!A.IsZero(k, i) || k == i)
        {
          ALU_ids.insert(std::make_pair(k, i));
          continue;
        }
        for (Index j = 0; j < i; ++j)
        {
          if (ALU_ids.contains(std::make_pair(k, j)) && ALU_ids.contains(std::make_pair(j, i)))
          {
            ALU_ids.insert(std::make_pair(k, i));
            break;
          }
        }
      }
    }
    auto ALU_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (const auto& pair : ALU_ids)
    {
      ALU_builder = ALU_builder.WithElement(pair.first, pair.second);
    }
    return SparseMatrixPolicy(ALU_builder, indexing_only);
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittleInPlace<SparseMatrixPolicy>::Decompose(SparseMatrixPolicy& ALU) const
  {
    const auto& views = views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrix::ViewType& alu_view) {
          auto aik_njk = views.aik_njk_.begin();
          auto aij_ajk = views.aij_ajk_.begin();
          auto aki_nji = views.aki_nji_.begin();
          auto akj_aji = views.akj_aji_.begin();
          for (const auto& nik_nki_aii : views.nik_nki_aii_)
          {
            for (Index ik = 0; ik < nik_nki_aii.first_; ++ik)
            {
              auto aik_view = alu_view.GetBlockView(aik_njk->first_);
              for (Index jk = 0; jk < aik_njk->second_; ++jk)
              {
                alu_view.ForEachBlock(
                    [](Real& aik, const Real& aij, const Real& ajk) { aik -= aij * ajk; },
                    aik_view,
                    alu_view.GetConstBlockView(aij_ajk->first_),
                    alu_view.GetConstBlockView(aij_ajk->second_));
                ++aij_ajk;
              }
              ++aik_njk;
            }
            for (Index ki = 0; ki < nik_nki_aii.second_; ++ki)
            {
              auto aki_view = alu_view.GetBlockView(aki_nji->first_);
              for (Index ji = 0; ji < aki_nji->second_; ++ji)
              {
                alu_view.ForEachBlock(
                    [](Real& aki, const Real& akj, const Real& aji) { aki -= akj * aji; },
                    aki_view,
                    alu_view.GetConstBlockView(akj_aji->first_),
                    alu_view.GetConstBlockView(akj_aji->second_));
                ++akj_aji;
              }
              alu_view.ForEachBlock(
                  [](Real& aki, const Real& aii) { aki /= aii; }, aki_view, alu_view.GetConstBlockView(nik_nki_aii.third_));
              ++aki_nji;
            }
          }
        },
        ALU)(ALU);
  }
}  // namespace micm
