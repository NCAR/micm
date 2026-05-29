// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  inline LuDecompositionMozartInPlace::LuDecompositionMozartInPlace()
  {
  }

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
    std::size_t n = matrix.NumRows();
    auto alu = GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value, true);
    for (std::size_t i = 0; i < n; ++i)
    {
      if (alu.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      std::tuple<std::size_t, std::size_t, std::size_t> aii_nji_nki(alu.VectorIndex(0, i, i), 0, 0);
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (alu.IsZero(j, i))
        {
          continue;
        }
        aji_.push_back(alu.VectorIndex(0, j, i));
        ++(std::get<1>(aii_nji_nki));
      }
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (alu.IsZero(i, k))
        {
          continue;
        }
        std::pair<std::size_t, std::size_t> aik_njk(alu.VectorIndex(0, i, k), 0);
        for (std::size_t j = i + 1; j < n; ++j)
        {
          if (alu.IsZero(j, i))
          {
            continue;
          }
          std::pair<std::size_t, std::size_t> ajk_aji(alu.VectorIndex(0, j, k), alu.VectorIndex(0, j, i));
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
    auto alu_builder = SparseMatrixPolicy::Create(n).SetNumberOfBlocks(A.NumberOfBlocks()).InitialValue(initial_value);
    for (auto& pair : ALU_ids)
    {
      alu_builder = alu_builder.WithElement(pair.first, pair.second);
    }
    return SparseMatrixPolicy(alu_builder, indexing_only);
  }

  template<class SparseMatrixPolicy>
    requires(!VectorizableSparse<SparseMatrixPolicy>)
  inline void LuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    const std::size_t N = ALU.NumRows();

    // Loop over blocks
    for (std::size_t i_block = 0; i_block < ALU.NumberOfBlocks(); ++i_block)
    {
      auto alu_vector = std::next(ALU.AsVector().begin(), i_block * ALU.FlatBlockSize());
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();

      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        const typename SparseMatrixPolicy::value_type Aii_inverse = 1.0 / alu_vector[std::get<0>(aii_nji_nki)];
        for (std::size_t ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          alu_vector[*aji] *= Aii_inverse;
          ++aji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const typename SparseMatrixPolicy::value_type Aik = alu_vector[std::get<0>(*aik_njk)];
          for (std::size_t ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
          {
            alu_vector[ajk_aji->first] -= alu_vector[ajk_aji->second] * Aik;
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
    const std::size_t N = ALU.NumRows();
    const std::size_t ALU_BLOCK_SIZE = ALU.NumberOfBlocks();
    constexpr std::size_t ALU_GROUP_VECTOR_SIZE = SparseMatrixPolicy::GroupVectorSize();
    const std::size_t ALU_GROUP_SIZE_OF_FLAT_BLOCK_SIZE = ALU.GroupSize();
    std::vector<double> Aii_inverse(ALU_GROUP_VECTOR_SIZE);

    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < ALU.NumberOfGroups(ALU_BLOCK_SIZE); ++i_group)
    {
      auto alu_vector = std::next(ALU.AsVector().begin(), i_group * ALU_GROUP_SIZE_OF_FLAT_BLOCK_SIZE);
      const std::size_t N_CELLS = std::min(ALU_GROUP_VECTOR_SIZE, ALU_BLOCK_SIZE - i_group * ALU_GROUP_VECTOR_SIZE);
      auto aji = aji_.begin();
      auto aik_njk = aik_njk_.begin();
      auto ajk_aji = ajk_aji_.begin();
      for (const auto& aii_nji_nki : aii_nji_nki_)
      {
        auto Aii_inverse_it = Aii_inverse.begin();
        auto ALU_vector_it = alu_vector + std::get<0>(aii_nji_nki);
        for (std::size_t i = 0; i < N_CELLS; ++i)
          *(Aii_inverse_it++) = 1.0 / *(ALU_vector_it++);
        for (std::size_t ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
        {
          auto ALU_vector_it = alu_vector + *aji;
          auto Aii_inverse_it = Aii_inverse.begin();
          for (std::size_t i = 0; i < N_CELLS; ++i)
            *(ALU_vector_it++) *= *(Aii_inverse_it++);
          ++aji;
        }
        for (std::size_t ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
        {
          const std::size_t aik = std::get<0>(*aik_njk);
          for (std::size_t ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
          {
            auto ALU_vector_first_it = alu_vector + ajk_aji->first;
            auto ALU_vector_second_it = alu_vector + ajk_aji->second;
            auto ALU_vector_aik_it = alu_vector + aik;
            for (std::size_t i = 0; i < N_CELLS; ++i)
              *(ALU_vector_first_it++) -= *(ALU_vector_second_it++) * *(ALU_vector_aik_it++);
            ++ajk_aji;
          }
          ++aik_njk;
        }
      }
    }
  }
}  // namespace micm
