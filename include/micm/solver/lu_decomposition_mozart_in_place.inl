// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  inline LuDecompositionMozartInPlace::LuDecompositionMozartInPlace() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace::LuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionMozartInPlace LuDecompositionMozartInPlace::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozartInPlace lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionMozartInPlace::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    Index n = matrix.NumRows();
    auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value, true);
    for (Index i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      std::tuple<Index, Index, Index> aii_nji_nki(ALU.VectorIndex(0, i, i), 0, 0);
      for (Index j = i + 1; j < n; ++j)
      {
        if (ALU.IsZero(j, i))
        {
          continue;
        }
        aji_.push_back(ALU.VectorIndex(0, j, i));
        ++(std::get<1>(aii_nji_nki));
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
        {
          continue;
        }
        std::pair<Index, Index> aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (Index j = i + 1; j < n; ++j)
        {
          if (ALU.IsZero(j, i))
          {
            continue;
          }
          std::pair<Index, Index> ajk_aji(ALU.VectorIndex(0, j, k), ALU.VectorIndex(0, j, i));
          ajk_aji_.push_back(ajk_aji);
          ++(std::get<1>(aik_njk));
        }
        aik_njk_.push_back(aik_njk);
        ++(std::get<2>(aii_nji_nki));
      }
      aii_nji_nki_.push_back(aii_nji_nki);
    }
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline SparseMatrixPolicy LuDecompositionMozartInPlace::GetLUMatrix(
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
  inline void LuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    SparseMatrixPolicy::Function(
        [this](auto&& alu_view)
        {
          auto aji = aji_.begin();
          auto aik_njk = aik_njk_.begin();
          auto ajk_aji = ajk_aji_.begin();
          auto Aii_inverse = alu_view.GetBlockVariable();
          for (const auto& aii_nji_nki : aii_nji_nki_)
          {
            alu_view.ForEachBlock(
                [](Real& aii_inv, const Real& aii) { aii_inv = 1.0 / aii; },
                Aii_inverse,
                alu_view.GetConstBlockView(std::get<0>(aii_nji_nki)));
            for (Index ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
            {
              alu_view.ForEachBlock(
                  [](Real& aji, const Real& aii_inv) { aji *= aii_inv; }, alu_view.GetBlockView(*aji), Aii_inverse);
              ++aji;
            }
            for (Index ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
            {
              auto aik_view = alu_view.GetBlockView(std::get<0>(*aik_njk));
              for (Index ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
              {
                alu_view.ForEachBlock(
                    [](Real& ajk, const Real& aji, const Real& aik) { ajk -= aji * aik; },
                    alu_view.GetBlockView(ajk_aji->first),
                    alu_view.GetConstBlockView(ajk_aji->second),
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
