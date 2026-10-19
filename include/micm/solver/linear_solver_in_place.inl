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
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy { return LuDecompositionPolicy::Create(m); })
  {
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolverInPlace(
      LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>&& other) noexcept
      : nLij_(std::move(other.nLij_)),
        Lij_yj_(std::move(other.Lij_yj_)),
        nUij_Uii_(std::move(other.nUij_Uii_)),
        Uij_xj_(std::move(other.Uij_xj_)),
        views_(nLij_, Lij_yj_, nUij_Uii_, Uij_xj_),
        lu_decomp_(std::move(other.lu_decomp_))
  {
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>&
  LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::operator=(
      LinearSolverInPlace<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>&& other) noexcept
  {
    if (this != &other)
    {
      nLij_ = std::move(other.nLij_);
      Lij_yj_ = std::move(other.Lij_yj_);
      nUij_Uii_ = std::move(other.nUij_Uii_);
      Uij_xj_ = std::move(other.Uij_xj_);
      views_ = Views(nLij_, Lij_yj_, nUij_Uii_, Uij_xj_);
      lu_decomp_ = std::move(other.lu_decomp_);
    }
    return *this;
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
        views_(),
        lu_decomp_(create_lu_decomp(matrix))
  {
    auto lu = lu_decomp_.GetLUMatrix(matrix, initial_value, true);
    std::vector<Index> nLij_temp;
    std::vector<IndexPair> Lij_yj_temp;
    std::vector<IndexPair> nUij_Uii_temp;
    std::vector<IndexPair> Uij_xj_temp;
    for (Index i = 0; i < lu.NumRows(); ++i)
    {
      Index nLij = 0;
      for (Index j = 0; j < i; ++j)
      {
        if (lu.IsZero(i, j))
        {
          continue;
        }
        Lij_yj_temp.push_back({ lu.VectorIndex(0, i, j), j });
        ++nLij;
      }
      nLij_temp.push_back(nLij);
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
        Uij_xj_temp.push_back({ lu.VectorIndex(0, i, j), j });
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_temp.push_back({ nUij, lu.VectorIndex(0, i, i) });
    }
    nLij_ = nLij_temp;
    Lij_yj_ = Lij_yj_temp;
    nUij_Uii_ = nUij_Uii_temp;
    Uij_xj_ = Uij_xj_temp;
    nLij_.CopyToDevice();
    Lij_yj_.CopyToDevice();
    nUij_Uii_.CopyToDevice();
    Uij_xj_.CopyToDevice();
    views_ = Views(nLij_, Lij_yj_, nUij_Uii_, Uij_xj_);
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
    const auto& views = views_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename DenseMatrix::ViewType& x_view, const typename SparseMatrix::ConstViewType& lu_view)
        {
          // For small vector lengths, the sums go into a temporary of L values that can stay in registers. The x
          // columns are in the same matrix, so the compiler cannot keep x_i in registers while it reads the other
          // x_j columns. On aarch64 with GCC 13, the temporary made the solve 1.4 times faster for L = 4, but 2 times
          // slower for L = 8 and L = 128, so for L > 4 the sums go directly into x_i.
          constexpr bool sum_in_registers = GROUP_VECTOR_SIZE_V<MatrixPolicy> <= 4;
          auto sum = lu_view.GetBlockVariable();

          // Subtracts count terms a_ij * x_j from x_i, with the iterator over the (a_ij, j) pairs
          auto subtract_terms = [&](auto& x_col_i, auto& a_x, const Index count)
          {
            if constexpr (sum_in_registers)
            {
              lu_view.ForEachBlock([](Real& xi_sum, const Real& xi) { xi_sum = xi; }, sum, x_col_i);
              for (Index k = 0; k < count; ++k)
              {
                lu_view.ForEachBlock(
                    [](Real& xi_sum, const Real& aij, const Real& xj) { xi_sum -= aij * xj; },
                    sum,
                    lu_view.GetConstBlockView((*a_x).first_),
                    x_view.GetConstColumnView((*a_x).second_));
                ++a_x;
              }
              lu_view.ForEachBlock([](Real& xi, const Real& xi_sum) { xi = xi_sum; }, x_col_i, sum);
            }
            else
            {
              for (Index k = 0; k < count; ++k)
              {
                lu_view.ForEachBlock(
                    [](Real& xi, const Real& aij, const Real& xj) { xi -= aij * xj; },
                    x_col_i,
                    lu_view.GetConstBlockView((*a_x).first_),
                    x_view.GetConstColumnView((*a_x).second_));
                ++a_x;
              }
            }
          };

          // Forward substitution
          // b values passed in as x; overwrites b with y values
          {
            auto Lij_yj = views.Lij_yj_.begin();
            Index i = 0;
            for (const auto& nLij : views.nLij_)
            {
              if (nLij > 0)
              {
                auto x_col_i = x_view.GetColumnView(i);
                subtract_terms(x_col_i, Lij_yj, nLij);
              }
              ++i;
            }
          }
          // Backward substitution
          // overwrites y values with x values
          {
            auto Uij_xj = views.Uij_xj_.begin();
            Index i = views.nUij_Uii_.size();
            for (const auto& nUij_Uii : views.nUij_Uii_)
            {
              --i;
              auto x_col_i = x_view.GetColumnView(i);
              subtract_terms(x_col_i, Uij_xj, nUij_Uii.first_);
              lu_view.ForEachBlock(
                  [](Real& xi, const Real& Uii) { xi /= Uii; }, x_col_i, lu_view.GetConstBlockView(nUij_Uii.second_));
            }
          }
        },
        x,
        lu_matrix)(x, lu_matrix);
}
}  // namespace micm