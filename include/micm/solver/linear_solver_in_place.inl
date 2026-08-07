// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

namespace micm
{
  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolverInPlace(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value)
      : LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy
            { return LuDecompositionPolicy::Create(m); })
  {
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolverInPlace(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value,
      const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)>& create_lu_decomp)
      : nLij_(),
        Lij_yj_(),
        nUij_Uii_(),
        Uij_xj_(),
        lu_decomp_(create_lu_decomp(matrix))
  {
    auto lu = lu_decomp_.GetLUMatrix(matrix, initial_value, true);
    for (Index i = 0; i < lu.NumRows(); ++i)
    {
      Index nLij = 0;
      for (Index j = 0; j < i; ++j)
      {
        if (lu.IsZero(i, j))
        {
          continue;
        }
        Lij_yj_.push_back(std::make_pair(lu.VectorIndex(0, i, j), j));
        ++nLij;
      }
      nLij_.push_back(nLij);
    }
    for (Index i = lu.NumRows() - 1; i != static_cast<Index>(-1); --i)
    {
      Index nUij = 0;
      for (Index j = i + 1; j < lu.NumColumns(); ++j)
      {
        if (lu.IsZero(i, j))
        {
          continue;
        }
        Uij_xj_.push_back(std::make_pair(lu.VectorIndex(0, i, j), j));
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_.push_back(std::make_pair(nUij, lu.VectorIndex(0, i, i)));
    }
  };

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      SparseMatrixPolicy& matrix) const
  {
    lu_decomp_.Decompose(matrix);
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lu_matrix) const
  {
    SparseMatrixPolicy::Function(
        [this](auto&& x_view, auto&& lu_view)
        {
          // Forward substitution
          // b values passed in as x; overwrites b with y values
          {
            auto Lij_yj = Lij_yj_.begin();
            Index i = 0;
            for (const auto& nLij : nLij_)
            {
              auto x_col_i = x_view.GetColumnView(i);
              for (Index k = 0; k < nLij; ++k)
              {
                lu_view.ForEachBlock(
                    [](Real& yi, const Real& Lij, const Real& yj) { yi -= Lij * yj; },
                    x_col_i,
                    lu_view.GetConstBlockView((*Lij_yj).first),
                    x_view.GetConstColumnView((*Lij_yj).second));
                ++Lij_yj;
              }
              ++i;
            }
          }
          // Backward substitution
          // overwrites y values with x values
          {
            auto Uij_xj = Uij_xj_.begin();
            Index i = nUij_Uii_.size();
            for (const auto& nUij_Uii : nUij_Uii_)
            {
              --i;
              auto x_col_i = x_view.GetColumnView(i);
              for (Index k = 0; k < nUij_Uii.first; ++k)
              {
                lu_view.ForEachBlock(
                    [](Real& xi, const Real& Uij, const Real& xj) { xi -= Uij * xj; },
                    x_col_i,
                    lu_view.GetConstBlockView((*Uij_xj).first),
                    x_view.GetConstColumnView((*Uij_xj).second));
                ++Uij_xj;
              }
              lu_view.ForEachBlock(
                  [](Real& xi, const Real& Uii) { xi /= Uii; }, x_col_i, lu_view.GetConstBlockView(nUij_Uii.second));
            }
          }
        },
        x,
        lu_matrix)(x, lu_matrix);
  }
}  // namespace micm