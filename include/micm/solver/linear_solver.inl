// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

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
  inline std::vector<std::size_t> DiagonalMarkowitzReorder(const MatrixPolicy& matrix)
  {
    const std::size_t order = matrix.NumRows();
    assert(order == matrix.NumColumns() && "Markowitz reorder requires a square matrix");
    // output_neighbors[v] = { c : edge v->c }, incoming_neighbors[c] = { v : edge v->c }, over the remaining nodes.
    std::vector<std::set<std::size_t>> output_neighbors(order), incoming_neighbors(order);
    for (std::size_t i = 0; i < order; ++i)
    {
      for (std::size_t j = 0; j < order; ++j)
      {
        if (matrix[i][j] != 0)
        {
          output_neighbors[i].insert(j);
          incoming_neighbors[j].insert(i);
        }
      }
    }
    std::vector<std::size_t> row_deg(order), col_deg(order);
    for (std::size_t v = 0; v < order; ++v)
    {
      row_deg[v] = output_neighbors[v].size();
      col_deg[v] = incoming_neighbors[v].size();
    }
    std::vector<char> alive(order, 1);
    std::vector<std::size_t> perm;
    perm.reserve(order);
    for (std::size_t step = 0; step < order; ++step)
    {
      // Select the remaining node with minimum Markowitz cost (row_deg-1)*(col_deg-1).
      // The diagonal keeps every live node's degrees >= 1, so the subtraction never underflows.
      std::size_t pivot = order;
      std::size_t best_cost = std::numeric_limits<std::size_t>::max();
      for (std::size_t v = 0; v < order; ++v)
      {
        if (!alive[v])
        {
          continue;
        }
        const std::size_t cost = (row_deg[v] - 1) * (col_deg[v] - 1);
        if (pivot == order || cost < best_cost)
        {
          best_cost = cost;
          pivot = v;
        }
      }
      perm.push_back(pivot);
      alive[pivot] = 0;
      std::vector<std::size_t> cols, ins;
      for (std::size_t c : output_neighbors[pivot])
      {
        if (c != pivot && alive[c])
        {
          cols.push_back(c);
        }
      }
      for (std::size_t i : incoming_neighbors[pivot])
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
      for (std::size_t i : ins)
      {
        for (std::size_t c : cols)
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
      for (std::size_t c : cols)
      {
        incoming_neighbors[c].erase(pivot);
        --col_deg[c];
      }
      for (std::size_t i : ins)
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
            { return LuDecompositionPolicy::template Create<SparseMatrixPolicy>(m); })
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
    auto lu = lu_decomp_.template GetLUMatrices<SparseMatrixPolicy>(matrix, initial_value, true);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    for (std::size_t i = 0; i < lower_matrix.NumRows(); ++i)
    {
      std::size_t nLij = 0;
      for (std::size_t j = 0; j < i; ++j)
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
    for (std::size_t i = upper_matrix.NumRows() - 1; i != static_cast<std::size_t>(-1); --i)
    {
      std::size_t nUij = 0;
      for (std::size_t j = i + 1; j < upper_matrix.NumColumns(); ++j)
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

    /// create and save the solver function
    auto x = MatrixPolicy(matrix.NumberOfBlocks(), matrix.NumRows());
    solve_func_ = LinearSolveFunc(x, lower_matrix, upper_matrix);
  };

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      const SparseMatrixPolicy& matrix,
      SparseMatrixPolicy& lower_matrix,
      SparseMatrixPolicy& upper_matrix) const
  {
    lu_decomp_.template Decompose<SparseMatrixPolicy>(matrix, lower_matrix, upper_matrix);
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lower_matrix,
      const SparseMatrixPolicy& upper_matrix) const
  {
    solve_func_(x, lower_matrix, upper_matrix);
  }

  template<class MatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline std::function<void(MatrixPolicy&, const SparseMatrixPolicy&, const SparseMatrixPolicy&)>
  LinearSolver<MatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolveFunc(
      const MatrixPolicy& x, const SparseMatrixPolicy& lower_matrix, const SparseMatrixPolicy& upper_matrix)
  {
    return SparseMatrixPolicy::Function(
      [nLij_Lii = nLij_Lii_,
       Lij_yj = Lij_yj_,
       nUij_Uii = nUij_Uii_,
       Uij_xj = Uij_xj_](auto&& x, auto&& lower, auto&& upper){
        // Forward Substitution
        // b values passed in as x; overwrites b values with y values
        {
          auto Lij_yj_it = Lij_yj.begin();
          std::size_t i = 0;
          for (const auto& nLij_Lii_entry : nLij_Lii)
          {
            for (std::size_t k = 0; k < nLij_Lii_entry.first; ++k)
            {
              lower.ForEachBlock(
                [&](double& yi, const double& Lij, const double& yj)
                {
                  yi -= Lij * yj;
                },
                x.GetColumnView(i),
                lower.GetConstBlockView((*Lij_yj_it).first),
                x.GetConstColumnView((*Lij_yj_it).second)
              );
              ++Lij_yj_it;
            }
            lower.ForEachBlock(
              [&](double& yi, const double& Lii)
              {
                yi /= Lii;
              },
              x.GetColumnView(i),
              lower.GetConstBlockView(nLij_Lii_entry.second)
            );
            ++i;
          }
        }
        // Backward Substitution
        // overwrites y values with x values
        {
          auto Uij_xj_it = Uij_xj.begin();
          std::size_t i = size(nUij_Uii);
          for (const auto& nUij_Uii_entry : nUij_Uii)
          {
            --i;
            for (std::size_t k = 0; k < nUij_Uii_entry.first; ++k)
            {
              upper.ForEachBlock(
                [&](double& xi, const double& Uij, const double& xj)
                {
                  xi -= Uij * xj;
                },
                x.GetColumnView(i),
                upper.GetConstBlockView((*Uij_xj_it).first),
                x.GetConstColumnView((*Uij_xj_it).second)
              );
              ++Uij_xj_it;
            }
            upper.ForEachBlock(
              [&](double& xi, const double& Uii)
              {
                xi /= Uii;
              },
              x.GetColumnView(i),
              upper.GetConstBlockView(nUij_Uii_entry.second)
            );
          }
        }
      },
      x,
      lower_matrix,
      upper_matrix);
  }
}  // namespace micm