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

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value)
      : LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy
            { return LuDecompositionPolicy::template Create<SparseMatrixPolicy>(m); })
  {
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
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
  };

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      const SparseMatrixPolicy& matrix,
      SparseMatrixPolicy& lower_matrix,
      SparseMatrixPolicy& upper_matrix) const
  {
    lu_decomp_.template Decompose<SparseMatrixPolicy>(matrix, lower_matrix, upper_matrix);
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<class MatrixPolicy>
    requires(!VectorizableDense<MatrixPolicy> || !VectorizableSparse<SparseMatrixPolicy>)
  inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lower_matrix,
      const SparseMatrixPolicy& upper_matrix) const
  {
    for (std::size_t i_cell = 0; i_cell < x.NumRows(); ++i_cell)
    {
      auto x_cell = x[i_cell];
      const std::size_t lower_grid_offset = i_cell * lower_matrix.FlatBlockSize();
      const std::size_t upper_grid_offset = i_cell * upper_matrix.FlatBlockSize();
      auto& y_cell = x_cell;  // Alias x for consistency with equations, but to reuse memory

      // Forward Substitution
      {
        auto y_elem = y_cell.begin();
        auto Lij_yj = Lij_yj_.begin();
        for (const auto& nLij_Lii : nLij_Lii_)
        {
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            *y_elem -= lower_matrix.AsVector()[lower_grid_offset + (*Lij_yj).first] * y_cell[(*Lij_yj).second];
            ++Lij_yj;
          }
          *(y_elem++) /= lower_matrix.AsVector()[lower_grid_offset + nLij_Lii.second];
        }
      }

      // Backward Substitution
      {
        auto x_elem = std::next(x_cell.end(), -1);
        auto Uij_xj = Uij_xj_.begin();
        for (const auto& nUij_Uii : nUij_Uii_)
        {
          // x_elem starts out as y_elem from the previous loop
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            *x_elem -= upper_matrix.AsVector()[upper_grid_offset + (*Uij_xj).first] * x_cell[(*Uij_xj).second];
            ++Uij_xj;
          }

          *(x_elem) /= upper_matrix.AsVector()[upper_grid_offset + nUij_Uii.second];
          // don't iterate before the beginning of the vector
          if (x_elem != x_cell.begin())
          {
            --x_elem;
          }
        }
      }
    }
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lower_matrix,
      const SparseMatrixPolicy& upper_matrix) const
  {
    constexpr std::size_t n_cells = MatrixPolicy::GroupVectorSize();
    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < x.NumberOfGroups(); ++i_group)
    {
      auto x_group = std::next(x.AsVector().begin(), i_group * x.GroupSize());
      auto L_group = std::next(lower_matrix.AsVector().begin(), i_group * lower_matrix.GroupSize());
      auto U_group = std::next(upper_matrix.AsVector().begin(), i_group * upper_matrix.GroupSize());
      // Forward Substitution
      {
        auto y_elem = x_group;
        auto Lij_yj = Lij_yj_.begin();
        for (const auto& nLij_Lii : nLij_Lii_)
        {
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            const std::size_t Lij_yj_first = (*Lij_yj).first;
            const std::size_t Lij_yj_second_times_n_cells = (*Lij_yj).second * n_cells;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              y_elem[i_cell] -= L_group[Lij_yj_first + i_cell] * x_group[Lij_yj_second_times_n_cells + i_cell];
            }
            ++Lij_yj;
          }
          const std::size_t nLij_Lii_second = nLij_Lii.second;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            y_elem[i_cell] /= L_group[nLij_Lii_second + i_cell];
          }
          y_elem += n_cells;
        }
      }

      // Backward Substitution
      {
        auto x_elem = std::next(x_group, x.GroupSize() - n_cells);
        auto Uij_xj = Uij_xj_.begin();
        for (const auto& nUij_Uii : nUij_Uii_)
        {
          // x_elem starts out as y_elem from the previous loop
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            const std::size_t Uij_xj_first = (*Uij_xj).first;
            const std::size_t Uij_xj_second_times_n_cells = (*Uij_xj).second * n_cells;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            {
              x_elem[i_cell] -= U_group[Uij_xj_first + i_cell] * x_group[Uij_xj_second_times_n_cells + i_cell];
            }
            ++Uij_xj;
          }
          const std::size_t nUij_Uii_second = nUij_Uii.second;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
          {
            x_elem[i_cell] /= U_group[nUij_Uii_second + i_cell];
          }

          // don't iterate before the beginning of the vector
          const std::size_t x_elem_distance = std::distance(x.AsVector().begin(), x_elem);
          x_elem -= std::min(n_cells, x_elem_distance);
        }
      }
    }
  }
}  // namespace micm