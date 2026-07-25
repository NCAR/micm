// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

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
    std::size_t n = matrix.NumRows();
    auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value, true);
    for (std::size_t i = 0; i < n; ++i)
    {
      if (ALU.IsZero(i, i))
      {
        throw std::runtime_error("Diagonal element is zero in LU decomposition");
      }
      std::tuple<std::size_t, std::size_t, std::size_t> aii_nji_nki(ALU.VectorIndex(0, i, i), 0, 0);
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (ALU.IsZero(j, i))
        {
          continue;
        }
        aji_.push_back(ALU.VectorIndex(0, j, i));
        ++(std::get<1>(aii_nji_nki));
      }
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (ALU.IsZero(i, k))
        {
          continue;
        }
        std::pair<std::size_t, std::size_t> aik_njk(ALU.VectorIndex(0, i, k), 0);
        for (std::size_t j = i + 1; j < n; ++j)
        {
          if (ALU.IsZero(j, i))
          {
            continue;
          }
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
      {
        if (!A.IsZero(i, j))
        {
          ALU_ids.insert(std::make_pair(i, j));
        }
      }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = i + 1; j < n; ++j)
      {
        if (ALU_ids.contains(std::make_pair(j, i)))
        {
          ALU_ids.insert(std::make_pair(j, i));
        }
      }
      for (std::size_t k = i + 1; k < n; ++k)
      {
        if (ALU_ids.contains(std::make_pair(i, k)))
        {
          for (std::size_t j = i + 1; j < n; ++j)
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
            [](double& aii_inv, const double& aii)
            {
              aii_inv = 1.0 / aii;
            },
            Aii_inverse,
            alu_view.GetConstBlockView(std::get<0>(aii_nji_nki)));
          for (std::size_t ij = 0; ij < std::get<1>(aii_nji_nki); ++ij)
          {
            alu_view.ForEachBlock(
              [](double& aji, const double& aii_inv)
              {
                aji *= aii_inv;
              },
              alu_view.GetBlockView(*aji),
              Aii_inverse);
            ++aji;
          }
          for (std::size_t ik = 0; ik < std::get<2>(aii_nji_nki); ++ik)
          {
            auto aik_view = alu_view.GetBlockView(std::get<0>(*aik_njk));
            for (std::size_t ijk = 0; ijk < std::get<1>(*aik_njk); ++ijk)
            {
              alu_view.ForEachBlock(
                [](double& ajk, const double& aji, const double& aik)
                {
                  ajk -= aji * aik;
                },
                alu_view.GetBlockView(ajk_aji->first),
                alu_view.GetConstBlockView(ajk_aji->second),
                aik_view);
              ++ajk_aji;
            }
            ++aik_njk;
          }
        }
      }
    )(ALU);
  } 
}  // namespace micm
