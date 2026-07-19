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
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittleInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    const Index n = ALU.NumRows();

    // Loop over blocks
    for (Index i_block = 0; i_block < ALU.NumberOfBlocks(); ++i_block)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_block * ALU.FlatBlockSize());
      auto aik_njk = aik_njk_.begin();
      auto aij_ajk = aij_ajk_.begin();
      auto aki_nji = aki_nji_.begin();
      auto akj_aji = akj_aji_.begin();

      for (const auto& nik_nki_aii : nik_nki_aii_)
      {
        for (Index ik = 0; ik < std::get<0>(nik_nki_aii); ++ik)
        {
          for (Index jk = 0; jk < aik_njk->second; ++jk)
          {
            ALU_vector[aik_njk->first] -= ALU_vector[aij_ajk->first] * ALU_vector[aij_ajk->second];
            ++aij_ajk;
          }
          ++aik_njk;
        }
        for (Index ki = 0; ki < std::get<1>(nik_nki_aii); ++ki)
        {
          for (Index ji = 0; ji < aki_nji->second; ++ji)
          {
            ALU_vector[aki_nji->first] -= ALU_vector[akj_aji->first] * ALU_vector[akj_aji->second];
            ++akj_aji;
          }
          ALU_vector[aki_nji->first] /= ALU_vector[std::get<2>(nik_nki_aii)];
          ++aki_nji;
        }
      }
    }
  }

  template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionDoolittleInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    const Index n = ALU.NumRows();
    const Index ALU_BlockSize = ALU.NumberOfBlocks();
    constexpr Index ALU_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const Index ALU_GroupSizeOfFlatBlockSize = ALU.GroupSize();

    // Loop over groups of blocks
    for (Index i_group = 0; i_group < ALU.NumberOfGroups(ALU_BlockSize); ++i_group)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_group * ALU_GroupSizeOfFlatBlockSize);
      const Index n_cells = std::min(ALU_GroupVectorSize, ALU_BlockSize - i_group * ALU_GroupVectorSize);
      auto aik_njk = aik_njk_.begin();
      auto aij_ajk = aij_ajk_.begin();
      auto aki_nji = aki_nji_.begin();
      auto akj_aji = akj_aji_.begin();
      for (const auto& nik_nki_aii : nik_nki_aii_)
      {
        const Index ik_limit = std::get<0>(nik_nki_aii);
        for (Index ik = 0; ik < ik_limit; ++ik)
        {
          const Index jk_limit = aik_njk->second;
          for (Index jk = 0; jk < jk_limit; ++jk)
          {
            auto ALU_vector_aik_njk_it = ALU_vector + aik_njk->first;
            auto ALU_vector_aij_ajk_first_it = ALU_vector + aij_ajk->first;
            auto ALU_vector_aij_ajk_second_it = ALU_vector + aij_ajk->second;
            for (Index i = 0; i < n_cells; ++i)
            {
              *(ALU_vector_aik_njk_it++) -= *(ALU_vector_aij_ajk_first_it++) * *(ALU_vector_aij_ajk_second_it++);
            }
            ++aij_ajk;
          }
          ++aik_njk;
        }
        const Index ki_limit = std::get<1>(nik_nki_aii);
        for (Index ki = 0; ki < ki_limit; ++ki)
        {
          const Index ji_limit = aki_nji->second;
          for (Index ji = 0; ji < ji_limit; ++ji)
          {
            auto ALU_vector_aki_nji_it = ALU_vector + aki_nji->first;
            auto ALU_vector_akj_aji_first_it = ALU_vector + akj_aji->first;
            auto ALU_vector_akj_aji_second_it = ALU_vector + akj_aji->second;
            for (Index i = 0; i < n_cells; ++i)
            {
              *(ALU_vector_aki_nji_it++) -= *(ALU_vector_akj_aji_first_it++) * *(ALU_vector_akj_aji_second_it++);
            }
            ++akj_aji;
          }
          auto ALU_vector_aki_nji_it = ALU_vector + aki_nji->first;
          auto ALU_vector_nik_nki_aii_it = ALU_vector + std::get<2>(nik_nki_aii);
          for (Index i = 0; i < n_cells; ++i)
          {
            *(ALU_vector_aki_nji_it++) /= *(ALU_vector_nik_nki_aii_it++);
          }
          ++aki_nji;
        }
      }
    }
  }

}  // namespace micm
