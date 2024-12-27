// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionMozartInPlace::LuDecompositionMozartInPlace()
  {
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline LuDecompositionMozartInPlace::LuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix)
  {
    Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline LuDecompositionMozartInPlace LuDecompositionMozartInPlace::Create(
      const SparseMatrixPolicy& matrix)
  {
    LuDecompositionMozartInPlace lu_decomp{};
    lu_decomp.Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());
    return lu_decomp;
  }

  template<class SparseMatrixPolicy>
  requires(SparseMatrixConcept<SparseMatrixPolicy>) inline void LuDecompositionMozartInPlace::Initialize(
      const SparseMatrixPolicy& matrix,
      auto initial_value)
  {
    MICM_PROFILE_FUNCTION();

    std::size_t n = matrix.NumRows();
    auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value);
    for (std::size_t i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i)) {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      std::tuple<std::size_t, std::size_t, std::size_t> aii_nji_nki(ALU.VectorIndex(0, i, i), 0, 0);
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (ALU.IsZero(j, i))
          continue;
        aji_.push_back(ALU.VectorIndex(0, j, i));
        ++(std::get<1>(aii_nji_nki));
      }
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
          continue;
        std::pair<std::size_t, std::size_t> aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (std::size_t j = i + 1; j < n; ++j)
        {
          if (ALU.IsZero(j, i))
            continue;
          std::pair<std::size_t, std::size_t> ajk_aji(ALU.VectorIndex(0, j, k), ALU.VectorIndex(0, j, i));
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
  requires(
      SparseMatrixConcept<SparseMatrixPolicy>) inline SparseMatrixPolicy LuDecompositionMozartInPlace::
      GetLUMatrix(const SparseMatrixPolicy& A, typename SparseMatrixPolicy::value_type initial_value)
  {
    MICM_PROFILE_FUNCTION();

    std::size_t n = A.NumRows();
    std::set<std::pair<std::size_t, std::size_t>> ALU_ids;
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j < n; ++j)
        if (!A.IsZero(i, j))
          ALU_ids.insert(std::make_pair(i, j));
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = i + 1; j < n; ++j)
        if (std::find(ALU_ids.begin(), ALU_ids.end(), std::make_pair(j, i)) != ALU_ids.end())
          ALU_ids.insert(std::make_pair(j, i));
      for (std::size_t k = i + 1; k < n; ++k)
        if (std::find(ALU_ids.begin(), ALU_ids.end(), std::make_pair(i, k)) != ALU_ids.end())
          for (std::size_t j = i + 1; j < n; ++j)
            if (std::find(ALU_ids.begin(), ALU_ids.end(), std::make_pair(j, i)) != ALU_ids.end())
              ALU_ids.insert(std::make_pair(j, k));
    }
    auto ALU_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : ALU_ids)
    {
      ALU_builder = ALU_builder.WithElement(pair.first, pair.second);
    }
    return SparseMatrixPolicy(ALU_builder);
  }

  template<class SparseMatrixPolicy>
  requires(!VectorizableSparse<SparseMatrixPolicy>) inline void LuDecompositionMozartInPlace::Decompose(
      SparseMatrixPolicy& ALU) const
  {
    MICM_PROFILE_FUNCTION();
    const std::size_t n = ALU.NumRows();

    // Loop over blocks
    for (std::size_t i_block = 0; i_block < ALU.NumberOfBlocks(); ++i_block)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_block * ALU.FlatBlockSize());
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();
      
      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        const typename SparseMatrixPolicy::value_type Aii_inverse = 1.0 / ALU_vector[std::get<0>(aii_nji_nki)];
        for (std::size_t ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          ALU_vector[*aji] *= Aii_inverse;
          ++aji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const typename SparseMatrixPolicy::value_type Aik = ALU_vector[std::get<0>(*aik_njk)];
          for (std::size_t ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
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
  requires(VectorizableSparse<SparseMatrixPolicy>) inline void LuDecompositionMozartInPlace::Decompose(
      SparseMatrixPolicy& ALU) const
  {
    MICM_PROFILE_FUNCTION();

    const std::size_t n = ALU.NumRows();
    const std::size_t ALU_BlockSize = ALU.NumberOfBlocks();
    constexpr std::size_t ALU_GroupVectorSize = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t ALU_GroupSizeOfFlatBlockSize = ALU.GroupSize();
    double Aii_inverse[ALU_GroupVectorSize];

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < ALU.NumberOfGroups(ALU_BlockSize); ++i_group)
    {
      auto ALU_vector = std::next(ALU.AsVector().begin(), i_group * ALU_GroupSizeOfFlatBlockSize);
      const std::size_t n_cells = std::min(ALU_GroupVectorSize, ALU_BlockSize - i_group * ALU_GroupVectorSize);
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();
      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        for (std::size_t i = 0; i < n_cells; ++i)
          Aii_inverse[i] = 1.0 / ALU_vector[std::get<0>(aii_nji_nki) + i];
        for (std::size_t ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          for (std::size_t i = 0; i < n_cells; ++i)
            ALU_vector[*aji + i] *= Aii_inverse[i];
          ++aji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const std::size_t aik = std::get<0>(*aik_njk);
          for (std::size_t ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
          {
            for (std::size_t i = 0; i < n_cells; ++i)
              ALU_vector[ajk_aji->first + i] -= ALU_vector[ajk_aji->second + i] * ALU_vector[aik + i];
            ++ajk_aji;
          }
          ++aik_njk;
        }
      }
    }
  }
}  // namespace micm
