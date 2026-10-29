// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

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
}  // namespace micm
