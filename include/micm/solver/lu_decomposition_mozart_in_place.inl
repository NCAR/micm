// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace<SparseMatrixPolicy>::LuDecompositionMozartInPlace() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace<SparseMatrixPolicy>::LuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix)
  {
    Initialize(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace<SparseMatrixPolicy>::LuDecompositionMozartInPlace(
      LuDecompositionMozartInPlace&& other) noexcept
      : aii_nji_nki_(std::move(other.aii_nji_nki_)),
        aji_(std::move(other.aji_)),
        aik_njk_(std::move(other.aik_njk_)),
        ajk_aji_(std::move(other.ajk_aji_)),
        views_(aii_nji_nki_, aji_, aik_njk_, ajk_aji_)
  {
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace<SparseMatrixPolicy>& LuDecompositionMozartInPlace<SparseMatrixPolicy>::operator=(
      LuDecompositionMozartInPlace&& other) noexcept
  {
    if (this != &other)
    {
      aii_nji_nki_ = std::move(other.aii_nji_nki_);
      aji_ = std::move(other.aji_);
      aik_njk_ = std::move(other.aik_njk_);
      ajk_aji_ = std::move(other.ajk_aji_);
      views_ = Views(aii_nji_nki_, aji_, aik_njk_, ajk_aji_);
    }
    return *this;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace<SparseMatrixPolicy> LuDecompositionMozartInPlace<SparseMatrixPolicy>::Create(
      const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozartInPlace<SparseMatrixPolicy> lu_decomp{};
    lu_decomp.Initialize(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozartInPlace<SparseMatrixPolicy>::Initialize(
      const SparseMatrixPolicy& matrix,
      auto initial_value)
  {
    Index n = matrix.NumRows();
    auto ALU = GetLUMatrix(matrix, initial_value, true);
    std::vector<IndexTrio> aii_nji_nki_temp;
    std::vector<Index> aji_temp;
    std::vector<IndexPair> aik_njk_temp;
    std::vector<IndexPair> ajk_aji_temp;
    for (Index i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      IndexTrio aii_nji_nki(ALU.VectorIndex(0, i, i), 0, 0);
      for (Index j = i + 1; j < n; ++j)
      {
        if (ALU.IsZero(j, i))
        {
          continue;
        }
        aji_temp.push_back(ALU.VectorIndex(0, j, i));
        ++(aii_nji_nki.second_);
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
        {
          continue;
        }
        IndexPair aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (Index j = i + 1; j < n; ++j)
        {
          if (ALU.IsZero(j, i))
          {
            continue;
          }
          IndexPair ajk_aji(ALU.VectorIndex(0, j, k), ALU.VectorIndex(0, j, i));
          ajk_aji_temp.push_back(ajk_aji);
          ++(aik_njk.second_);
        }
        aik_njk_temp.push_back(aik_njk);
        ++(aii_nji_nki.third_);
      }
      aii_nji_nki_temp.push_back(aii_nji_nki);
    }
    aii_nji_nki_ = aii_nji_nki_temp;
    aji_ = aji_temp;
    aik_njk_ = aik_njk_temp;
    ajk_aji_ = ajk_aji_temp;
    aii_nji_nki_.CopyToDevice();
    aji_.CopyToDevice();
    aik_njk_.CopyToDevice();
    ajk_aji_.CopyToDevice();
    views_ = Views(aii_nji_nki_, aji_, aik_njk_, ajk_aji_);
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline SparseMatrixPolicy LuDecompositionMozartInPlace<SparseMatrixPolicy>::GetLUMatrix(
      const SparseMatrixPolicy& A,
      typename SparseMatrixPolicy::value_type initial_value,
      bool indexing_only)
  {
    Index n = A.NumRows();
    std::set<std::pair<Index, Index>> ALU_ids;
    for (Index i = 0; i < n; ++i)
    {
      for (Index j = 0; j < n; ++j)
      {
        if (!A.IsZero(i, j))
        {
          ALU_ids.insert(std::make_pair(i, j));
        }
      }
    }
    for (Index i = 0; i < n; ++i)
    {
      for (Index j = i + 1; j < n; ++j)
      {
        if (ALU_ids.contains(std::make_pair(j, i)))
        {
          ALU_ids.insert(std::make_pair(j, i));
        }
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (ALU_ids.contains(std::make_pair(i, k)))
        {
          for (Index j = i + 1; j < n; ++j)
          {
            if (ALU_ids.contains(std::make_pair(j, i)))
            {
              ALU_ids.insert(std::make_pair(j, k));
            }
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
  inline void LuDecompositionMozartInPlace<SparseMatrixPolicy>::Decompose(SparseMatrixPolicy& ALU) const
  {
    const auto& views = views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(typename SparseMatrix::ViewType alu_view) {
          auto aji = views.aji_.begin();
          auto aik_njk = views.aik_njk_.begin();
          auto ajk_aji = views.ajk_aji_.begin();
          auto Aii_inverse = alu_view.GetBlockVariable();
          for (const auto& aii_nji_nki : views.aii_nji_nki_)
          {
            alu_view.ForEachBlock(
                [](Real& aii_inv, const Real& aii) { aii_inv = 1.0 / aii; },
                Aii_inverse,
                alu_view.GetConstBlockView(aii_nji_nki.first_));
            for (Index ij = 0; ij < aii_nji_nki.second_; ++ij)
            {
              alu_view.ForEachBlock(
                  [](Real& aji, const Real& aii_inv) { aji *= aii_inv; }, alu_view.GetBlockView(*aji), Aii_inverse);
              ++aji;
            }
            for (Index ik = 0; ik < aii_nji_nki.third_; ++ik)
            {
              auto aik_view = alu_view.GetBlockView((*aik_njk).first_);
              for (Index ijk = 0; ijk < (*aik_njk).second_; ++ijk)
              {
                alu_view.ForEachBlock(
                    [](Real& ajk, const Real& aji, const Real& aik) { ajk -= aji * aik; },
                    alu_view.GetBlockView(ajk_aji->first_),
                    alu_view.GetConstBlockView(ajk_aji->second_),
                    aik_view);
                ++ajk_aji;
              }
              ++aik_njk;
            }
          }
        },
        ALU)(ALU);
  }
}  // namespace micm
