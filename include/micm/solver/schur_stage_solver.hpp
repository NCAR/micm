// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition_doolittle.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <vector>

namespace micm
{

  /// @brief Schur-reduced stage solver for semi-explicit index-1 DAE Rosenbrock steps.
  ///
  /// The stage matrix of a mass-matrix Rosenbrock step is
  ///   [ alpha*I - J_xx    -J_xz ]   =   [ A_xx  A_xz ]
  ///   [     -J_zx         -J_zz ]       [ A_zx  A_zz ]
  /// with x the differential and z the algebraic variables. Index 1 guarantees
  /// A_zz nonsingular, so each stage solve reduces to the differential block:
  ///
  ///   1. per-group factor A_zz (the algebraic rows decompose into connected
  ///      components under the A_zz sparsity; each component is treated as a
  ///      small dense block),
  ///   2. form the Schur complement S = A_xx - A_xz * W with W = A_zz^{-1} A_zx
  ///      (numeric per step; the sparsity of S is computed once, symbolically,
  ///      as pattern(A_xx) + the product pattern through each group),
  ///   3. factor S with MICM's sparse Doolittle LU, and per stage solve
  ///        y_z = A_zz^{-1} b_z,   x = S^{-1} (b_x - A_xz y_z),   z = y_z - W x.
  ///
  /// The reduction is exact linear algebra: solutions match the unreduced
  /// [alpha*M - J] solve to roundoff. The factored dimension drops from
  /// n_x + n_z to n_x, which is the measured remedy for constraint-heavy
  /// systems (see docs/superpowers/notes/2026-07-20-schur-reduction-design.md).
  ///
  /// Symbolic setup happens once per instance from the Jacobian sparsity and
  /// the mass-matrix diagonal; Factor() runs the numeric phase per step
  /// attempt and returns false if any algebraic block is numerically singular
  /// (an index defect, e.g. a constraint row that does not involve its own
  /// algebraic variable). Host-side, standard-ordering matrices only.
  template<class SparseMatrixPolicy, class DenseMatrixPolicy>
  class SchurStageSolver
  {
   public:
    SchurStageSolver(const SparseMatrixPolicy& jacobian, const std::vector<double>& mass_matrix_diagonal)
        : n_cells_(jacobian.NumberOfBlocks()),
          jac_stride_(jacobian.FlatBlockSize())
    {
      const std::size_t n_vars = mass_matrix_diagonal.size();
      // Partition the state.
      std::vector<std::ptrdiff_t> x_slot_of(n_vars, -1), z_slot_of(n_vars, -1);
      for (std::size_t i = 0; i < n_vars; ++i)
      {
        if (mass_matrix_diagonal[i] != 0.0)
        {
          x_slot_of[i] = static_cast<std::ptrdiff_t>(x_of_.size());
          x_of_.push_back(i);
        }
        else
        {
          z_slot_of[i] = static_cast<std::ptrdiff_t>(z_of_.size());
          z_of_.push_back(i);
        }
      }
      n_x_ = x_of_.size();
      n_z_ = z_of_.size();

      // Enumerate the Jacobian pattern into the four blocks (local indices +
      // flat offsets into block 0 of the Jacobian data).
      struct RawEntry
      {
        std::size_t r, c, offset;
      };
      std::vector<RawEntry> zz_raw;
      std::vector<std::vector<RawEntry>> zx_raw_by_row(n_z_);
      const auto& row_start = jacobian.RowStart();
      const auto& row_ids = jacobian.RowIds();
      for (std::size_t r = 0; r < n_vars; ++r)
      {
        for (std::size_t k = row_start[r]; k < row_start[r + 1]; ++k)
        {
          const std::size_t c = row_ids[k];
          const std::size_t offset = k;  // CSR position == block-0 flat offset
          if (x_slot_of[r] >= 0 && x_slot_of[c] >= 0)
          {
            xx_entries_.push_back({ static_cast<std::size_t>(x_slot_of[r]), static_cast<std::size_t>(x_slot_of[c]), offset });
          }
          else if (x_slot_of[r] >= 0)
          {
            xz_entries_.push_back({ static_cast<std::size_t>(x_slot_of[r]), static_cast<std::size_t>(z_slot_of[c]), offset });
          }
          else if (x_slot_of[c] >= 0)
          {
            zx_raw_by_row[static_cast<std::size_t>(z_slot_of[r])].push_back(
                { static_cast<std::size_t>(z_slot_of[r]), static_cast<std::size_t>(x_slot_of[c]), offset });
          }
          else
          {
            zz_raw.push_back({ static_cast<std::size_t>(z_slot_of[r]), static_cast<std::size_t>(z_slot_of[c]), offset });
          }
        }
      }

      // Group algebraic rows into connected components of the A_zz pattern
      // (union-find), so each component factors as a small dense block.
      std::vector<std::size_t> parent(n_z_);
      std::iota(parent.begin(), parent.end(), 0);
      auto find = [&](std::size_t a)
      {
        while (parent[a] != a)
        {
          parent[a] = parent[parent[a]];
          a = parent[a];
        }
        return a;
      };
      for (const auto& e : zz_raw)
        parent[find(e.r)] = find(e.c);

      std::vector<std::ptrdiff_t> group_id(n_z_, -1);
      group_of_z_.resize(n_z_);
      z_pos_in_group_.resize(n_z_);
      for (std::size_t z = 0; z < n_z_; ++z)
      {
        const std::size_t root = find(z);
        if (group_id[root] < 0)
        {
          group_id[root] = static_cast<std::ptrdiff_t>(groups_.size());
          groups_.emplace_back();
        }
        Group& group = groups_[static_cast<std::size_t>(group_id[root])];
        group_of_z_[z] = static_cast<std::size_t>(group_id[root]);
        z_pos_in_group_[z] = group.z_slots.size();
        group.z_slots.push_back(z);
      }
      for (const auto& e : zz_raw)
      {
        Group& group = groups_[group_of_z_[e.r]];
        group.zz_entries.push_back({ z_pos_in_group_[e.r], z_pos_in_group_[e.c], e.offset });
      }
      for (std::size_t z = 0; z < n_z_; ++z)
      {
        Group& group = groups_[group_of_z_[z]];
        for (const auto& e : zx_raw_by_row[z])
        {
          std::size_t col_index = group.x_cols.size();
          for (std::size_t j = 0; j < group.x_cols.size(); ++j)
          {
            if (group.x_cols[j] == e.c)
            {
              col_index = j;
              break;
            }
          }
          if (col_index == group.x_cols.size())
            group.x_cols.push_back(e.c);
          group.zx_entries.push_back({ z_pos_in_group_[z], col_index, e.offset });
        }
      }
      std::size_t dense_size = 0, w_size = 0, pivot_base = 0;
      for (auto& group : groups_)
      {
        group.dense_offset = dense_size;
        group.w_offset = w_size;
        pivot_base_.push_back(pivot_base);
        dense_size += group.z_slots.size() * group.z_slots.size();
        w_size += group.z_slots.size() * group.x_cols.size();
        pivot_base += group.z_slots.size();
      }
      zz_lu_.resize(n_cells_ * dense_size);
      zz_dense_stride_ = dense_size;
      zz_piv_.resize(n_cells_ * n_z_);
      w_.resize(n_cells_ * w_size);
      w_stride_ = w_size;
      xz_values_.resize(n_cells_ * xz_entries_.size());
      yz_.resize(n_z_);

      // Symbolic Schur complement: pattern(A_xx) + diagonal + the product
      // pattern of A_xz through each group's x-column set.
      auto builder = SparseMatrixPolicy::Create(n_x_).SetNumberOfBlocks(n_cells_).InitialValue(0.0);
      for (std::size_t r = 0; r < n_x_; ++r)
        builder = builder.WithElement(r, r);
      for (const auto& e : xx_entries_)
        builder = builder.WithElement(e.r, e.c);
      for (const auto& e : xz_entries_)
      {
        const Group& group = groups_[group_of_z_[e.c]];
        for (const std::size_t col : group.x_cols)
          builder = builder.WithElement(e.r, col);
      }
      s_ = SparseMatrixPolicy(builder);
      s_stride_ = s_.FlatBlockSize();

      // Precompute flat scatter offsets into S (block 0).
      for (std::size_t r = 0; r < n_x_; ++r)
        s_diag_offsets_.push_back(s_.VectorIndex(0, r, r));
      for (const auto& e : xx_entries_)
        s_xx_offsets_.push_back(s_.VectorIndex(0, e.r, e.c));
      for (std::size_t i = 0; i < xz_entries_.size(); ++i)
      {
        const auto& e = xz_entries_[i];
        const Group& group = groups_[group_of_z_[e.c]];
        const std::size_t zpos = z_pos_in_group_[e.c];
        for (std::size_t j = 0; j < group.x_cols.size(); ++j)
        {
          products_.push_back({ i,
                                group.w_offset + zpos * group.x_cols.size() + j,
                                s_.VectorIndex(0, e.r, group.x_cols[j]) });
        }
      }

      // Sparse LU machinery for S.
      auto lu_pair = LuDecompositionDoolittle::GetLUMatrices<SparseMatrixPolicy>(s_, 0.0);
      s_lower_ = std::move(lu_pair.first);
      s_upper_ = std::move(lu_pair.second);
      s_solver_ = LinearSolver<SparseMatrixPolicy, LuDecompositionDoolittle>(s_, 0.0);
      x_rhs_ = DenseMatrixPolicy(n_cells_, n_x_, 0.0);
    }

    std::size_t NumDifferential() const
    {
      return n_x_;
    }

    std::size_t NumAlgebraic() const
    {
      return n_z_;
    }

    /// @brief Numeric phase for one step attempt. `negative_jacobian` holds -J
    ///        (unshifted); `alpha` = 1/(h*gamma). Returns false if an algebraic
    ///        block is numerically singular.
    bool Factor(const SparseMatrixPolicy& negative_jacobian, double alpha)
    {
      const auto& jac = negative_jacobian.AsVector();
      auto& s_data = s_.AsVector();

      for (std::size_t cell = 0; cell < n_cells_; ++cell)
      {
        const std::size_t jac_base = cell * jac_stride_;
        // Cache A_xz values for this cell (reused by S formation and stage solves).
        for (std::size_t i = 0; i < xz_entries_.size(); ++i)
          xz_values_[cell * xz_entries_.size() + i] = jac[jac_base + xz_entries_[i].offset];

        // Per-group dense factorization of A_zz and W = A_zz^{-1} A_zx.
        for (std::size_t g = 0; g < groups_.size(); ++g)
        {
          const Group& group = groups_[g];
          const std::size_t nz = group.z_slots.size();
          double* a = zz_lu_.data() + cell * zz_dense_stride_ + group.dense_offset;
          int* piv = zz_piv_.data() + cell * n_z_ + pivot_base_[g];
          std::fill(a, a + nz * nz, 0.0);
          for (const auto& e : group.zz_entries)
            a[e.r * nz + e.c] = jac[jac_base + e.offset];
          if (!DenseLuFactor(a, piv, nz))
            return false;

          double* w = w_.data() + cell * w_stride_ + group.w_offset;
          const std::size_t ncols = group.x_cols.size();
          std::fill(w, w + nz * ncols, 0.0);
          for (const auto& e : group.zx_entries)
            w[e.r * ncols + e.c] = jac[jac_base + e.offset];
          // Solve A_zz W = A_zx column-wise (W stored row-major nz x ncols).
          DenseLuSolveMultiple(a, piv, nz, w, ncols);
        }

        // Numeric Schur complement: S = alpha*I + (-J)_xx - A_xz * W.
        const std::size_t s_base = cell * s_stride_;
        std::fill(s_data.begin() + s_base, s_data.begin() + s_base + s_stride_, 0.0);
        for (std::size_t r = 0; r < n_x_; ++r)
          s_data[s_base + s_diag_offsets_[r]] += alpha;
        for (std::size_t i = 0; i < xx_entries_.size(); ++i)
          s_data[s_base + s_xx_offsets_[i]] += jac[jac_base + xx_entries_[i].offset];
        const double* xz_vals = xz_values_.data() + cell * xz_entries_.size();
        const double* w_cell = w_.data() + cell * w_stride_;
        for (const auto& op : products_)
          s_data[s_base + op.s_offset] -= xz_vals[op.xz_index] * w_cell[op.w_index];
      }

      s_solver_.Factor(s_, s_lower_, s_upper_);
      return true;
    }

    /// @brief In-place stage solve for all cells: b holds the stage right-hand
    ///        side on entry and the stage solution on exit.
    void Solve(DenseMatrixPolicy& b)
    {
      // Eliminate z from the right-hand side.
      for (std::size_t cell = 0; cell < n_cells_; ++cell)
      {
        for (std::size_t g = 0; g < groups_.size(); ++g)
        {
          const Group& group = groups_[g];
          const std::size_t nz = group.z_slots.size();
          double* yz = yz_.data() + pivot_base_[g];
          for (std::size_t i = 0; i < nz; ++i)
            yz[i] = b[cell][z_of_[group.z_slots[i]]];
          const double* a = zz_lu_.data() + cell * zz_dense_stride_ + group.dense_offset;
          DenseLuSolveMultiple(a, zz_piv_.data() + cell * n_z_ + pivot_base_[g], nz, yz, 1);
        }
        for (std::size_t r = 0; r < n_x_; ++r)
          x_rhs_[cell][r] = b[cell][x_of_[r]];
        const double* xz_vals = xz_values_.data() + cell * xz_entries_.size();
        for (std::size_t i = 0; i < xz_entries_.size(); ++i)
        {
          const std::size_t z = xz_entries_[i].c;
          x_rhs_[cell][xz_entries_[i].r] -= xz_vals[i] * yz_[pivot_base_[group_of_z_[z]] + z_pos_in_group_[z]];
        }
        // Stash y_z for the recovery pass (b's z slots are free to reuse).
        for (std::size_t z = 0; z < n_z_; ++z)
          b[cell][z_of_[z]] = yz_[pivot_base_[group_of_z_[z]] + z_pos_in_group_[z]];
      }

      // Differential solve on the Schur complement (all cells at once).
      s_solver_.Solve(x_rhs_, s_lower_, s_upper_);

      // Recover z = y_z - W x and scatter the solution back.
      for (std::size_t cell = 0; cell < n_cells_; ++cell)
      {
        const double* w_cell = w_.data() + cell * w_stride_;
        for (const auto& group : groups_)
        {
          const std::size_t nz = group.z_slots.size();
          const std::size_t ncols = group.x_cols.size();
          const double* w = w_cell + group.w_offset;
          for (std::size_t i = 0; i < nz; ++i)
          {
            double correction = 0.0;
            for (std::size_t j = 0; j < ncols; ++j)
              correction += w[i * ncols + j] * x_rhs_[cell][group.x_cols[j]];
            b[cell][z_of_[group.z_slots[i]]] -= correction;
          }
        }
        for (std::size_t r = 0; r < n_x_; ++r)
          b[cell][x_of_[r]] = x_rhs_[cell][r];
      }
    }

   private:
    struct BlockEntry
    {
      std::size_t r, c, offset;  // local row, local col, flat Jacobian offset (block 0)
    };
    struct Group
    {
      std::vector<std::size_t> z_slots;      // z slots in this connected component
      std::vector<std::size_t> x_cols;       // union of coupled differential columns
      std::vector<BlockEntry> zz_entries;    // r, c are group-local positions
      std::vector<BlockEntry> zx_entries;    // r group-local, c indexes x_cols
      std::size_t dense_offset = 0;          // into zz_lu_ per cell
      std::size_t w_offset = 0;              // into w_ per cell
    };
    struct ProductOp
    {
      std::size_t xz_index;  // into xz_entries_/xz_values_
      std::size_t w_index;   // into w_ per cell
      std::size_t s_offset;  // flat offset into S (block 0)
    };

    /// Dense LU with partial pivoting; returns false on a (near-)zero pivot.
    static bool DenseLuFactor(double* a, int* piv, std::size_t n)
    {
      for (std::size_t k = 0; k < n; ++k)
      {
        std::size_t p = k;
        double biggest = std::abs(a[k * n + k]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
          if (std::abs(a[i * n + k]) > biggest)
          {
            biggest = std::abs(a[i * n + k]);
            p = i;
          }
        }
        if (!(biggest > 0.0) || !std::isfinite(biggest))
          return false;
        piv[k] = static_cast<int>(p);
        if (p != k)
        {
          for (std::size_t j = 0; j < n; ++j)
            std::swap(a[k * n + j], a[p * n + j]);
        }
        const double inv = 1.0 / a[k * n + k];
        for (std::size_t i = k + 1; i < n; ++i)
        {
          const double factor = a[i * n + k] * inv;
          a[i * n + k] = factor;
          for (std::size_t j = k + 1; j < n; ++j)
            a[i * n + j] -= factor * a[k * n + j];
        }
      }
      return true;
    }

    /// Solve A X = B for row-major X (n x m) using the factors from DenseLuFactor.
    static void DenseLuSolveMultiple(const double* a, const int* piv, std::size_t n, double* x, std::size_t m)
    {
      for (std::size_t k = 0; k < n; ++k)
      {
        const std::size_t p = static_cast<std::size_t>(piv[k]);
        if (p != k)
        {
          for (std::size_t j = 0; j < m; ++j)
            std::swap(x[k * m + j], x[p * m + j]);
        }
        for (std::size_t i = k + 1; i < n; ++i)
        {
          const double factor = a[i * n + k];
          for (std::size_t j = 0; j < m; ++j)
            x[i * m + j] -= factor * x[k * m + j];
        }
      }
      for (std::size_t k = n; k-- > 0;)
      {
        const double inv = 1.0 / a[k * n + k];
        for (std::size_t j = 0; j < m; ++j)
        {
          double sum = x[k * m + j];
          for (std::size_t i = k + 1; i < n; ++i)
            sum -= a[k * n + i] * x[i * m + j];
          x[k * m + j] = sum * inv;
        }
      }
    }

    std::size_t n_cells_;
    std::size_t jac_stride_;
    std::size_t n_x_ = 0, n_z_ = 0;
    std::vector<std::size_t> x_of_, z_of_;
    std::vector<std::size_t> group_of_z_, z_pos_in_group_;
    std::vector<std::size_t> pivot_base_;
    std::vector<BlockEntry> xx_entries_, xz_entries_;
    std::vector<Group> groups_;
    std::vector<ProductOp> products_;
    std::vector<std::size_t> s_diag_offsets_, s_xx_offsets_;
    SparseMatrixPolicy s_;
    std::size_t s_stride_ = 0, zz_dense_stride_ = 0, w_stride_ = 0;
    SparseMatrixPolicy s_lower_, s_upper_;
    LinearSolver<SparseMatrixPolicy, LuDecompositionDoolittle> s_solver_;
    DenseMatrixPolicy x_rhs_;
    std::vector<double> zz_lu_, w_, xz_values_, yz_;
    std::vector<int> zz_piv_;
  };

}  // namespace micm
