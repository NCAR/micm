// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{

  inline LuDecompositionDoolittleInPlace::LuDecompositionDoolittleInPlace() = default;

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace::LuDecompositionDoolittleInPlace(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline LuDecompositionDoolittleInPlace LuDecompositionDoolittleInPlace::Create(const SparseMatrixPolicy& matrix)
  {
    LuDecompositionDoolittleInPlace lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittleInPlace::Initialize(const SparseMatrixPolicy& matrix, auto initial_value)
  {
    Index n = matrix.NumRows();
    auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value, true);
    for (Index i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      std::tuple<Index, Index, Index> nik_nki_aii(0, 0, ALU.VectorIndex(0, i, i));
      for (Index k = i; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
        {
          continue;
        }
        std::pair<Index, Index> aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (Index j = 0; j < i; ++j)
        {
          if (ALU.IsZero(i, j) || ALU.IsZero(j, k))
          {
            continue;
          }
          aij_ajk_.push_back(std::make_pair(ALU.VectorIndex(0, i, j), ALU.VectorIndex(0, j, k)));
          ++(aik_njk.second);
        }
        aik_njk_.push_back(aik_njk);
        ++(std::get<0>(nik_nki_aii));
      }
      for (Index k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(k, i))
        {
          continue;
        }
        std::pair<Index, Index> aki_nji(ALU.VectorIndex(0, k, i), 0);
        for (Index j = 0; j < i; ++j)
        {
          if (ALU.IsZero(k, j) || ALU.IsZero(j, i))
          {
            continue;
          }
          akj_aji_.push_back(std::make_pair(ALU.VectorIndex(0, k, j), ALU.VectorIndex(0, j, i)));
          ++(aki_nji.second);
        }
        aki_nji_.push_back(aki_nji);
        ++(std::get<1>(nik_nki_aii));
      }
      nik_nki_aii_.push_back(nik_nki_aii);
    }
  }

  template<class SparseMatrixPolicy>
    requires(SparseMatrixConcept<SparseMatrixPolicy>)
  inline SparseMatrixPolicy LuDecompositionDoolittleInPlace::GetLUMatrix(
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
  inline void LuDecompositionDoolittleInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    SparseMatrixPolicy::Function(
        [this](auto&& alu_view)
        {
          auto aik_njk = aik_njk_.begin();
          auto aij_ajk = aij_ajk_.begin();
          auto aki_nji = aki_nji_.begin();
          auto akj_aji = akj_aji_.begin();
          for (const auto& nik_nki_aii : nik_nki_aii_)
          {
            for (Index ik = 0; ik < std::get<0>(nik_nki_aii); ++ik)
            {
              auto aik_view = alu_view.GetBlockView(aik_njk->first);
              for (Index jk = 0; jk < aik_njk->second; ++jk)
              {
                alu_view.ForEachBlock(
                    [](Real& aik, const Real& aij, const Real& ajk) { aik -= aij * ajk; },
                    aik_view,
                    alu_view.GetConstBlockView(aij_ajk->first),
                    alu_view.GetConstBlockView(aij_ajk->second));
                ++aij_ajk;
              }
              ++aik_njk;
            }
            for (Index ki = 0; ki < std::get<1>(nik_nki_aii); ++ki)
            {
              auto aki_view = alu_view.GetBlockView(aki_nji->first);
              for (Index ji = 0; ji < aki_nji->second; ++ji)
              {
                alu_view.ForEachBlock(
                    [](Real& aki, const Real& akj, const Real& aji) { aki -= akj * aji; },
                    aki_view,
                    alu_view.GetConstBlockView(akj_aji->first),
                    alu_view.GetConstBlockView(akj_aji->second));
                ++akj_aji;
              }
              alu_view.ForEachBlock(
                  [](Real& aki, const Real& aii) { aki /= aii; },
                  aki_view,
                  alu_view.GetConstBlockView(std::get<2>(nik_nki_aii)));
              ++aki_nji;
            }
          }
        },
        ALU)(ALU);
  }
}  // namespace micm
