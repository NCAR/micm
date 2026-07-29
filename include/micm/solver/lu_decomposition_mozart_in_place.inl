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
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    const Index n = ALU.NumRows();

    // Loop over blocks
    for (Index i_block = 0; i_block < ALU.NumberOfBlocks(); ++i_block)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_block * ALU.FlatBlockSize());
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();

      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        const typename SparseMatrixPolicy::value_type Aii_inverse = 1.0 / ALU_vector[std::get<0>(aii_nji_nki)];
        for (Index ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          ALU_vector[*aji] *= Aii_inverse;
          ++aji;
        }
        for (Index ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const typename SparseMatrixPolicy::value_type Aik = ALU_vector[std::get<0>(*aik_njk)];
          for (Index ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
          {
            ALU_vector[ajk_aji->first] -= ALU_vector[ajk_aji->second] * Aik;
            ++ajk_aji;
          }
          ++aik_njk;
        }
      }
    }
  }

  template<class SparseMatrixPolicy>
    requires(VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    const Index n = ALU.NumRows();
    const Index ALU_BlockSize = ALU.NumberOfBlocks();
    constexpr Index ALU_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const Index ALU_GroupSizeOfFlatBlockSize = ALU.GroupSize();
    std::vector<Real> Aii_inverse(ALU_GroupVectorSize);

    // Loop over groups of blocks
    for (Index i_group = 0; i_group < ALU.NumberOfGroups(ALU_BlockSize); ++i_group)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_group * ALU_GroupSizeOfFlatBlockSize);
      const Index n_cells = std::min(ALU_GroupVectorSize, ALU_BlockSize - i_group * ALU_GroupVectorSize);
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();
      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        auto Aii_inverse_it = Aii_inverse.begin();
        auto ALU_vector_it = ALU_vector + std::get<0>(aii_nji_nki);
        for (Index i = 0; i < n_cells; ++i)
        {
          *(Aii_inverse_it++) = 1.0 / *(ALU_vector_it++);
        }
        for (Index ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          auto ALU_vector_it = ALU_vector + *aji;
          auto Aii_inverse_it = Aii_inverse.begin();
          for (Index i = 0; i < n_cells; ++i)
          {
            *(ALU_vector_it++) *= *(Aii_inverse_it++);
          }
          ++aji;
        }
        for (Index ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const Index aik = std::get<0>(*aik_njk);
          for (Index ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
          {
            auto ALU_vector_first_it = ALU_vector + ajk_aji->first;
            auto ALU_vector_second_it = ALU_vector + ajk_aji->second;
            auto ALU_vector_aik_it = ALU_vector + aik;
            for (Index i = 0; i < n_cells; ++i)
            {
              *(ALU_vector_first_it++) -= *(ALU_vector_second_it++) * *(ALU_vector_aik_it++);
            }
            ++ajk_aji;
          }
          ++aik_njk;
        }
      }
    }
  }
}  // namespace micm
