// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/types.hpp>

#include <limits>
#include <set>
#include <vector>

namespace micm
{
  // Diagonal Markowitz (minimum-degree) reordering on a SPARSE adjacency representation.
  //
  // This is the same algorithm as the classic dense Markowitz scan -- it selects, at each
  // elimination step, the remaining pivot with the lowest Markowitz cost (row_deg-1)*(col_deg-1)
  // and accounts for the fill that eliminating it introduces -- and it produces an ordering of
  // identical minimum-degree quality. Only the representation differs, which is why it looks
  // nothing like the dense version:
  //   * The dense version mutates an order x order pattern matrix, recomputes every candidate's
  //     row/column degree by re-scanning the trailing submatrix each step, and physically swaps
  //     rows/columns to move the chosen pivot into place -- O(order^3) time, O(order^2) memory.
  //   * This version stores the sparsity pattern as a directed graph (edge v->c iff
  //     matrix[v][c] != 0), keeps the degrees up to date incrementally, and marks eliminated
  //     nodes with an `alive` flag instead of moving anything -- ~O(order^2 + fill) time, with
  //     memory proportional to nonzeros + fill. This is essential for large mechanisms (thousands
  //     of species), where the dense form is intractable.
  //
  // Returns perm where perm[new_index] = old_index.
  template<class MatrixPolicy>
  inline std::vector<Index> DiagonalMarkowitzReorder(const MatrixPolicy& matrix)
  {
    const Index order = matrix.NumRows();
    assert(order == matrix.NumColumns() && "Markowitz reorder requires a square matrix");
    // output_neighbors[v] = { c : edge v->c }, incoming_neighbors[c] = { v : edge v->c }, over the remaining nodes.
    std::vector<std::set<Index>> output_neighbors(order), incoming_neighbors(order);
    for (Index i = 0; i < order; ++i)
    {
      for (Index j = 0; j < order; ++j)
      {
        if (matrix[i][j] != 0)
        {
          output_neighbors[i].insert(j);
          incoming_neighbors[j].insert(i);
        }
      }
    }
    std::vector<Index> row_deg(order), col_deg(order);
    for (Index v = 0; v < order; ++v)
    {
      row_deg[v] = output_neighbors[v].size();
      col_deg[v] = incoming_neighbors[v].size();
    }
    std::vector<Bool> alive(order, 1);
    std::vector<Index> perm;
    perm.reserve(order);
    for (Index step = 0; step < order; ++step)
    {
      // Select the remaining node with minimum Markowitz cost (row_deg-1)*(col_deg-1).
      // The diagonal keeps every live node's degrees >= 1, so the subtraction never underflows.
      Index pivot = order;
      Index best_cost = std::numeric_limits<Index>::max();
      for (Index v = 0; v < order; ++v)
      {
        if (!alive[v])
        {
          continue;
        }
        const Index cost = (row_deg[v] - 1) * (col_deg[v] - 1);
        if (pivot == order || cost < best_cost)
        {
          best_cost = cost;
          pivot = v;
        }
      }
      perm.push_back(pivot);
      alive[pivot] = 0;
      std::vector<Index> cols, ins;
      for (Index c : output_neighbors[pivot])
      {
        if (c != pivot && alive[c])
        {
          cols.push_back(c);
        }
      }
      for (Index i : incoming_neighbors[pivot])
      {
        if (i != pivot && alive[i])
        {
          ins.push_back(i);
        }
      }
      // Fill: eliminating the pivot couples everything that pointed into it with everything it
      // pointed out to. Concretely, for every live in-neighbor i (edge i->pivot) and every live
      // out-neighbor c (edge pivot->c), a new edge i->c must exist in the factored matrix. This
      // is the sparse equivalent of the dense version's "OR the pivot row down each column"
      // fill step. insert(c).second is true only when the edge is genuinely new, so degrees are
      // bumped exactly once per introduced fill element -- keeping row_deg/col_deg exact without
      // any rescan.
      for (Index i : ins)
      {
        for (Index c : cols)
        {
          if (output_neighbors[i].insert(c).second)
          {
            ++row_deg[i];
            incoming_neighbors[c].insert(i);
            ++col_deg[c];
          }
        }
      }
      // Drop the eliminated pivot from its live neighbors' degree counts.
      for (Index c : cols)
      {
        incoming_neighbors[c].erase(pivot);
        --col_deg[c];
      }
      for (Index i : ins)
      {
        output_neighbors[i].erase(pivot);
        --row_deg[i];
      }
    }
    return perm;
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value)
      : LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy
            { return LuDecompositionPolicy::Create(m); })
  {
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value,
      const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)>& create_lu_decomp)
      : nLij_Lii_(),
        Lij_yj_(),
        nUij_Uii_(),
        Uij_xj_(),
        lu_decomp_(create_lu_decomp(matrix))
  {
    auto lu = lu_decomp_.GetLUMatrices(matrix, initial_value, true);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    for (Index i = 0; i < lower_matrix.NumRows(); ++i)
    {
      Index nLij = 0;
      for (Index j = 0; j < i; ++j)
      {
        if (lower_matrix.IsZero(i, j))
        {
          continue;
        }
        Lij_yj_.push_back(std::make_pair(lower_matrix.VectorIndex(0, i, j), j));
        ++nLij;
      }
      // There must always be a non-zero element on the diagonal
      nLij_Lii_.push_back(std::make_pair(nLij, lower_matrix.VectorIndex(0, i, i)));
    }
    for (Index i = upper_matrix.NumRows() - 1; i != static_cast<Index>(-1); --i)
    {
      Index nUij = 0;
      for (Index j = i + 1; j < upper_matrix.NumColumns(); ++j)
      {
        if (upper_matrix.IsZero(i, j))
        {
          continue;
        }
        Uij_xj_.push_back(std::make_pair(upper_matrix.VectorIndex(0, i, j), j));
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_.push_back(std::make_pair(nUij, upper_matrix.VectorIndex(0, i, i)));
    }
  };

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      const SparseMatrixPolicy& matrix,
      SparseMatrixPolicy& lower_matrix,
      SparseMatrixPolicy& upper_matrix) const
  {
    lu_decomp_.Decompose(matrix, lower_matrix, upper_matrix);
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lower_matrix,
      const SparseMatrixPolicy& upper_matrix) const
  {
    SparseMatrixPolicy::Function(
        [this](auto&& x_view, auto&& lower_view, auto&& upper_view)
        {
          // Forward Substitution
          // b values passed in as x; overwrites b values with y values
          {
            auto Lij_yj = Lij_yj_.begin();
            Index i = 0;
            for (const auto& nLij_Lii : nLij_Lii_)
            {
              auto x_col_i = x_view.GetColumnView(i);
              for (Index k = 0; k < nLij_Lii.first; ++k)
              {
                lower_view.ForEachBlock(
                    [](Real& yi, const Real& Lij, const Real& yj) { yi -= Lij * yj; },
                    x_col_i,
                    lower_view.GetConstBlockView((*Lij_yj).first),
                    x_view.GetConstColumnView((*Lij_yj).second));
                ++Lij_yj;
              }
              lower_view.ForEachBlock(
                  [](Real& yi, const Real& Lii) { yi /= Lii; }, x_col_i, lower_view.GetConstBlockView(nLij_Lii.second));
              ++i;
            }
          }
          // Backward Substitution
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
                upper_view.ForEachBlock(
                    [](Real& xi, const Real& Uij, const Real& xj) { xi -= Uij * xj; },
                    x_col_i,
                    upper_view.GetConstBlockView((*Uij_xj).first),
                    x_view.GetConstColumnView((*Uij_xj).second));
                ++Uij_xj;
              }
              upper_view.ForEachBlock(
                  [](Real& xi, const Real& Uii) { xi /= Uii; }, x_col_i, upper_view.GetConstBlockView(nUij_Uii.second));
            }
          }
        },
        x,
        lower_matrix,
        upper_matrix)(x, lower_matrix, upper_matrix);
  }
}  // namespace micm