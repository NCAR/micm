// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Unit tests for the "grouped" view enrichment introduced to eliminate per-
// element index recomputation inside the ForEachBlock / ForEachRow hot path.
//
// The grouped views (`GroupedColumnView` and `GroupedConstColumnView` on the
// dense matrices, `GroupedBlockView` and `GroupedConstBlockView` on the sparse
// orderings) carry a precomputed base_ pointer scoped to a specific group. They
// are returned by `GetColumnView` / `GetConstColumnView` / `GetBlockView` /
// `GetConstBlockView` when those methods are called on a `GroupView` /
// `ConstGroupView` (as opposed to the raw matrix).
//
// These tests cover:
//   1. That the grouped views carry the correct view-category tag.
//   2. That element access (read + write) via ForEachBlock / ForEachRow yields
//      the same values whether we pass raw views or grouped views.
//   3. That the fast-path GetBlockElement / GetRowElement overloads for the
//      grouped tag are dispatched to (verified by relying on their behaviour
//      for a mixed dense + sparse call, which the fallback overloads compute
//      via different arithmetic paths).
//   4. That the L=1 degenerate case (plain Matrix, VectorMatrix<1>, standard-
//      ordered SparseMatrix) works identically to L>1.
//   5. That the const variant of a mutable GroupView also produces a grouped
//      view (via GetConstColumnView / GetConstBlockView on a mutable GroupView).

#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering_compressed_sparse_column.hpp>
#include <micm/util/sparse_matrix_standard_ordering_compressed_sparse_row.hpp>
#include <micm/util/sparse_matrix_vector_ordering_compressed_sparse_column.hpp>
#include <micm/util/sparse_matrix_vector_ordering_compressed_sparse_row.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/view_category.hpp>

#include <gtest/gtest.h>

#include <type_traits>

using namespace micm;

// ---------------------------------------------------------------------------
// Category tag checks (compile-time; no runtime work).
// ---------------------------------------------------------------------------

TEST(GroupedView, CategoryTagsAreDistinct)
{
  // The four category tags are separate types; conflating any of them would
  // route element accesses to the wrong GetBlockElement overload.
  static_assert(!std::is_same_v<GroupedDenseMatrixColumnViewTag, DenseMatrixColumnViewTag>);
  static_assert(!std::is_same_v<GroupedSparseMatrixBlockViewTag, SparseMatrixBlockViewTag>);
  static_assert(!std::is_same_v<GroupedDenseMatrixColumnViewTag, GroupedSparseMatrixBlockViewTag>);

  // The grouped concepts must not accidentally match the raw view concepts
  // and vice-versa; if they did, both overloads would be viable and dispatch
  // would be ambiguous.
  using GCV = Matrix<double>::GroupView::GroupedColumnView;
  static_assert(GroupedDenseMatrixColumnView<GCV>);
  static_assert(!DenseMatrixColumnView<GCV>);
  static_assert(!SparseMatrixBlockView<GCV>);
  static_assert(!GroupedSparseMatrixBlockView<GCV>);
  static_assert(!VectorLike<GCV>);

  using GBV = SparseMatrix<double>::GroupView::GroupedBlockView;
  static_assert(GroupedSparseMatrixBlockView<GBV>);
  static_assert(!SparseMatrixBlockView<GBV>);
  static_assert(!DenseMatrixColumnView<GBV>);
  static_assert(!GroupedDenseMatrixColumnView<GBV>);
  static_assert(!VectorLike<GBV>);
}

// ---------------------------------------------------------------------------
// Dense matrix (L=1): GroupView.GetColumnView returns a grouped view whose
// base_ pointer maps to the exact element and round-trips through ForEachRow.
// ---------------------------------------------------------------------------

TEST(GroupedView, MatrixColumnViewRoundTrip)
{
  Matrix<double> matrix{ 3, 4, 0.0 };
  // Set distinct values so aliasing bugs would show up immediately.
  for (std::size_t r = 0; r < 3; ++r)
  {
    for (std::size_t c = 0; c < 4; ++c)
    {
      matrix[r][c] = static_cast<double>(10 * r + c);
    }
  }

  // Row 1 of the matrix acts as the "group".
  Matrix<double>::GroupView row_view(matrix, 1);
  auto col2_mut = row_view.GetColumnView(2);
  auto col2_const = row_view.GetConstColumnView(2);
  auto col3_const = row_view.GetConstColumnView(3);

  static_assert(GroupedDenseMatrixColumnView<decltype(col2_mut)>);
  static_assert(GroupedDenseMatrixColumnView<decltype(col2_const)>);
  static_assert(GroupedDenseMatrixColumnView<decltype(col3_const)>);

  // Base pointer aims at the exact element (row=1, col=2, expected value 12).
  EXPECT_EQ(*col2_mut.base_, 12.0);
  EXPECT_EQ(*col2_const.base_, 12.0);
  EXPECT_EQ(*col3_const.base_, 13.0);

  // Write via ForEachRow: col2 = col3 * 2  ==> row 1 only.
  row_view.ForEachRow(
      [](double& a, const double& b) { a = b * 2.0; }, col2_mut, col3_const);
  EXPECT_EQ(matrix[1][2], 26.0);
  EXPECT_EQ(matrix[0][2], 2.0);   // unchanged (still 2 from init)
  EXPECT_EQ(matrix[2][2], 22.0);  // unchanged (still 22 from init)
}

// ---------------------------------------------------------------------------
// Dense matrix (L>1): GroupView.GetColumnView base_ points at the L-row block.
// Element access is contiguous: base_[0..L-1].
// ---------------------------------------------------------------------------

TEST(GroupedView, VectorMatrixColumnViewContiguousBlock)
{
  // Two full L=4 groups of rows, 3 columns.
  VectorMatrix<double, 4> matrix{ 8, 3, 0.0 };
  for (std::size_t r = 0; r < 8; ++r)
  {
    for (std::size_t c = 0; c < 3; ++c)
    {
      matrix[r][c] = static_cast<double>(100 + 10 * r + c);
    }
  }

  // Group 1 spans rows 4..7.
  VectorMatrix<double, 4>::GroupView group1(matrix, 1);
  auto col1_mut = group1.GetColumnView(1);
  auto col2_const = group1.GetConstColumnView(2);

  static_assert(GroupedDenseMatrixColumnView<decltype(col1_mut)>);
  static_assert(GroupedDenseMatrixColumnView<decltype(col2_const)>);

  // base_ points at row 4 for the requested column; base_[i] is row 4+i.
  for (std::size_t i = 0; i < 4; ++i)
  {
    EXPECT_EQ(col1_mut.base_[i], static_cast<double>(100 + 10 * (4 + i) + 1));
    EXPECT_EQ(col2_const.base_[i], static_cast<double>(100 + 10 * (4 + i) + 2));
  }

  // ForEachRow scales col1 by 2 over the group; other groups untouched.
  group1.ForEachRow([](double& a, const double& b) { a = b * 2.0; }, col1_mut, col2_const);
  for (std::size_t r = 4; r < 8; ++r)
  {
    EXPECT_EQ(matrix[r][1], 2.0 * static_cast<double>(100 + 10 * r + 2));
  }
  // Rows outside the group are still what we initialised.
  for (std::size_t r = 0; r < 4; ++r)
  {
    EXPECT_EQ(matrix[r][1], static_cast<double>(100 + 10 * r + 1));
  }
}

// ---------------------------------------------------------------------------
// Dense matrix: raw and grouped views produce the same result when driven by
// ForEachRow through Function(). If the grouped path were broken (or the
// fallback path were being taken for a grouped view) the two would diverge.
// ---------------------------------------------------------------------------

template<template<class, std::size_t> class VMPolicy, std::size_t L>
void RunDenseGroupedVsRawEquivalence()
{
  using M = VMPolicy<double, L>;
  const std::size_t rows = 3 * L;
  const std::size_t cols = 4;
  M raw{ rows, cols, 0.0 };
  M grouped{ rows, cols, 0.0 };

  // Fill both matrices identically.
  for (std::size_t r = 0; r < rows; ++r)
  {
    for (std::size_t c = 0; c < cols; ++c)
    {
      raw[r][c] = static_cast<double>(r * 10 + c);
      grouped[r][c] = static_cast<double>(r * 10 + c);
    }
  }

  // "grouped" path: uses group view GetColumnView (returns grouped view).
  auto grouped_func = M::Function(
      [](auto&& m)
      {
        m.ForEachRow(
            [](const double& a, double& out) { out = a * 3.0 + 1.0; },
            m.GetConstColumnView(0),
            m.GetColumnView(3));
      },
      grouped);
  grouped_func(grouped);

  // "raw" path: same math, but read from the matrix directly. Any dispatch
  // bug in the grouped overloads would cause the two matrices to diverge.
  for (std::size_t r = 0; r < rows; ++r)
  {
    raw[r][3] = raw[r][0] * 3.0 + 1.0;
  }

  for (std::size_t r = 0; r < rows; ++r)
  {
    for (std::size_t c = 0; c < cols; ++c)
    {
      EXPECT_EQ(raw[r][c], grouped[r][c]) << "row=" << r << " col=" << c << " L=" << L;
    }
  }
}

TEST(GroupedView, DenseGroupedMatchesRawL1)
{
  RunDenseGroupedVsRawEquivalence<VectorMatrix, 1>();
}
TEST(GroupedView, DenseGroupedMatchesRawL2)
{
  RunDenseGroupedVsRawEquivalence<VectorMatrix, 2>();
}
TEST(GroupedView, DenseGroupedMatchesRawL4)
{
  RunDenseGroupedVsRawEquivalence<VectorMatrix, 4>();
}
TEST(GroupedView, DenseGroupedMatchesRawL8)
{
  RunDenseGroupedVsRawEquivalence<VectorMatrix, 8>();
}

// Plain Matrix (L=1) as a separate template specialization path.
TEST(GroupedView, PlainMatrixGroupedRoundTrip)
{
  Matrix<double> raw{ 4, 3, 0.0 };
  Matrix<double> grouped{ 4, 3, 0.0 };
  for (std::size_t r = 0; r < 4; ++r)
  {
    for (std::size_t c = 0; c < 3; ++c)
    {
      raw[r][c] = static_cast<double>(r + c * 5);
      grouped[r][c] = static_cast<double>(r + c * 5);
    }
  }

  auto func = Matrix<double>::Function(
      [](auto&& m)
      {
        m.ForEachRow(
            [](const double& a, const double& b, double& out) { out = a - b; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            m.GetColumnView(2));
      },
      grouped);
  func(grouped);
  for (std::size_t r = 0; r < 4; ++r)
  {
    raw[r][2] = raw[r][0] - raw[r][1];
  }

  for (std::size_t r = 0; r < 4; ++r)
  {
    for (std::size_t c = 0; c < 3; ++c)
    {
      EXPECT_EQ(raw[r][c], grouped[r][c]) << "row=" << r << " col=" << c;
    }
  }
}

// ---------------------------------------------------------------------------
// Sparse matrix (standard, L=1): grouped block view carries a base_ + offset.
// ---------------------------------------------------------------------------

TEST(GroupedView, SparseStandardBlockViewRoundTrip)
{
  auto builder =
      SparseMatrix<double>::Create(3).WithElement(0, 1).WithElement(1, 2).WithElement(2, 2).SetNumberOfBlocks(2);
  SparseMatrix<double> matrix{ builder };
  matrix[0][0][1] = 1.0;
  matrix[0][1][2] = 2.0;
  matrix[0][2][2] = 3.0;
  matrix[1][0][1] = 10.0;
  matrix[1][1][2] = 20.0;
  matrix[1][2][2] = 30.0;

  // Group 1 = block 1 for standard ordering (L=1).
  SparseMatrix<double>::GroupView group1(matrix, 1);
  auto v01_mut = group1.GetBlockView(matrix.VectorIndex(0, 0, 1));
  auto v12_const = group1.GetConstBlockView(matrix.VectorIndex(0, 1, 2));
  static_assert(GroupedSparseMatrixBlockView<decltype(v01_mut)>);
  static_assert(GroupedSparseMatrixBlockView<decltype(v12_const)>);

  // group_base_ points at the start of block 1's data slice; offset picks the
  // right non-zero element. For standard ordering block_in_group is always 0.
  EXPECT_EQ(v01_mut.group_base_[v01_mut.block_offset_], 10.0);
  EXPECT_EQ(v12_const.group_base_[v12_const.block_offset_], 20.0);

  group1.ForEachBlock(
      [](double& a, const double& b) { a = b + 5.0; }, v01_mut, v12_const);
  EXPECT_EQ(matrix[1][0][1], 25.0);
  EXPECT_EQ(matrix[0][0][1], 1.0);  // untouched
}

// ---------------------------------------------------------------------------
// Sparse matrix (vector ordering, L=4): grouped block view spans L blocks
// contiguously. Element access is `group_base_[block_offset_ + block_in_group]`.
// ---------------------------------------------------------------------------

TEST(GroupedView, SparseVectorBlockViewContiguousBlock)
{
  using SM = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  auto builder = SM::Create(3).WithElement(0, 1).WithElement(1, 2).SetNumberOfBlocks(8);
  SM matrix{ builder };
  // Fill each block's non-zeros with a distinct pattern.
  for (std::size_t b = 0; b < 8; ++b)
  {
    matrix[b][0][1] = 100.0 + b;
    matrix[b][1][2] = 200.0 + b;
  }

  // Group 1 spans blocks 4..7 (L=4).
  SM::GroupView group1(matrix, 1);
  auto v01 = group1.GetBlockView(matrix.VectorIndex(0, 0, 1));
  auto v12 = group1.GetConstBlockView(matrix.VectorIndex(0, 1, 2));
  static_assert(GroupedSparseMatrixBlockView<decltype(v01)>);
  static_assert(GroupedSparseMatrixBlockView<decltype(v12)>);

  // group_base_ + block_offset_ + i lands on the (block_in_group = i) block.
  for (std::size_t i = 0; i < 4; ++i)
  {
    EXPECT_EQ(v01.group_base_[v01.block_offset_ + i], 100.0 + (4 + i));
    EXPECT_EQ(v12.group_base_[v12.block_offset_ + i], 200.0 + (4 + i));
  }

  // ForEachBlock touches only the L=4 blocks of the group.
  group1.ForEachBlock(
      [](double& a, const double& b) { a = a + b; }, v01, v12);
  for (std::size_t b = 4; b < 8; ++b)
  {
    EXPECT_EQ(matrix[b][0][1], (100.0 + b) + (200.0 + b));
  }
  for (std::size_t b = 0; b < 4; ++b)
  {
    EXPECT_EQ(matrix[b][0][1], 100.0 + b);  // other group untouched
  }
}

// ---------------------------------------------------------------------------
// Sparse + dense mixed inside a sparse GroupView's ForEachBlock: the sparse
// arg dispatches to the grouped-sparse GetBlockElement overload, the dense
// arg dispatches to the grouped-dense one. Verifies the two grouped tags
// coexist correctly for the linear-solver's actual call pattern.
// ---------------------------------------------------------------------------

// Forward declaration; body appears further down (needs the same
// SparseMatrix<...>::VectorIndex mapping the live matrix uses).
template<std::size_t L>
static std::size_t SparseRefIndexA();

template<std::size_t L>
void RunSparseVectorMixedGroupedEquivalence()
{
  using Sparse = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<L>>;
  using Dense = VectorMatrix<double, L>;

  const std::size_t blocks = 2 * L;
  auto builder = Sparse::Create(3).WithElement(0, 1).SetNumberOfBlocks(blocks);
  Sparse sparse_ref{ builder };
  Sparse sparse{ builder };
  Dense dense_ref{ blocks, 3, 0.0 };
  Dense dense{ blocks, 3, 0.0 };
  for (std::size_t b = 0; b < blocks; ++b)
  {
    sparse_ref[b][0][1] = 1.0 + b;
    sparse[b][0][1] = 1.0 + b;
    for (std::size_t c = 0; c < 3; ++c)
    {
      dense_ref[b][c] = 10.0 + b + c;
      dense[b][c] = 10.0 + b + c;
    }
  }

  // Same computation, "grouped" path via Function():
  //     dense[:, 0] = dense[:, 1] + sparse[:, 0, 1]
  auto func = Sparse::Function(
      [](auto&& sm, auto&& dm)
      {
        // sm is a sparse GroupView; dm is a dense GroupView.
        // GetConstBlockView / GetConstColumnView / GetColumnView on group views
        // return grouped views. The sparse GetBlockElement dispatches on the
        // *grouped* tags for both.
        sm.ForEachBlock(
            [](double& out, const double& a, const double& b) { out = a + b; },
            dm.GetColumnView(0),
            dm.GetConstColumnView(1),
            sm.GetConstBlockView(SparseRefIndexA<L>()));
      },
      sparse,
      dense);
  func(sparse, dense);

  // Reference: same math done by hand.
  for (std::size_t b = 0; b < blocks; ++b)
  {
    dense_ref[b][0] = dense_ref[b][1] + sparse_ref[b][0][1];
  }

  for (std::size_t b = 0; b < blocks; ++b)
  {
    for (std::size_t c = 0; c < 3; ++c)
    {
      EXPECT_EQ(dense[b][c], dense_ref[b][c]) << "b=" << b << " c=" << c << " L=" << L;
    }
  }
}

// Helper: fetch the sparse (0,1) vector index for the given L. We use a lambda
// to construct a throwaway matrix and read its VectorIndex mapping so we don't
// depend on internal ordering-policy details in the test.
template<std::size_t L>
static std::size_t SparseRefIndexA()
{
  using SM = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<L>>;
  static const std::size_t idx = []
  {
    auto b = SM::Create(3).WithElement(0, 1).SetNumberOfBlocks(2 * L);
    SM m{ b };
    return m.VectorIndex(0, 0, 1);
  }();
  return idx;
}

TEST(GroupedView, MixedSparseAndDenseGroupedL1)
{
  RunSparseVectorMixedGroupedEquivalence<1>();
}
TEST(GroupedView, MixedSparseAndDenseGroupedL2)
{
  RunSparseVectorMixedGroupedEquivalence<2>();
}
TEST(GroupedView, MixedSparseAndDenseGroupedL4)
{
  RunSparseVectorMixedGroupedEquivalence<4>();
}
TEST(GroupedView, MixedSparseAndDenseGroupedL8)
{
  RunSparseVectorMixedGroupedEquivalence<8>();
}

// ---------------------------------------------------------------------------
// The sparse GroupView also exposes (row, col) overloads for GetConstBlockView
// and GetBlockView. Those keep the raw ConstBlockView return type (they don't
// know the block_offset_ until VectorIndexFromRowColumn resolves the element),
// and remain routed through the SparseMatrixBlockView overload of
// GetBlockElement. This test guards against accidental breakage of that path
// while the grouped path is added alongside.
// ---------------------------------------------------------------------------

TEST(GroupedView, SparseVectorRowColOverloadStillFunctional)
{
  using SM = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  auto builder = SM::Create(3).WithElement(0, 1).SetNumberOfBlocks(8);
  SM matrix{ builder };
  for (std::size_t b = 0; b < 8; ++b)
  {
    matrix[b][0][1] = static_cast<double>(b + 1);
  }

  SM::ConstGroupView group1(matrix, 1);
  auto raw = group1.GetConstBlockView(0, 1);
  // The (row, col) overload returns the *raw* ConstBlockView from the parent
  // SparseMatrix, not the grouped one. Read via ForEachBlock and confirm.
  static_assert(SparseMatrixBlockView<decltype(raw)>);
  static_assert(!GroupedSparseMatrixBlockView<decltype(raw)>);

  double accumulator = 0.0;
  group1.ForEachBlock([&accumulator](const double& v) { accumulator += v; }, raw);
  // Group 1 covers blocks 4..7, values 5..8. Sum = 5+6+7+8 = 26.
  EXPECT_DOUBLE_EQ(accumulator, 26.0);
}

// ---------------------------------------------------------------------------
// Grouped views on a ConstGroupView also work (i.e. the const-only overloads
// exist and don't accidentally require a mutable parent group view).
// ---------------------------------------------------------------------------

TEST(GroupedView, ConstGroupViewProducesGroupedViews)
{
  using SM = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  auto builder = SM::Create(3).WithElement(0, 1).SetNumberOfBlocks(4);
  SM matrix{ builder };
  for (std::size_t b = 0; b < 4; ++b)
  {
    matrix[b][0][1] = static_cast<double>(b + 1);
  }

  const SM& cmatrix = matrix;
  SM::ConstGroupView cgv(cmatrix, 0);
  auto v = cgv.GetConstBlockView(matrix.VectorIndex(0, 0, 1));
  static_assert(GroupedSparseMatrixBlockView<decltype(v)>);
  for (std::size_t i = 0; i < 4; ++i)
  {
    EXPECT_EQ(v.group_base_[v.block_offset_ + i], static_cast<double>(i + 1));
  }

  VectorMatrix<double, 4> dense{ 4, 2, 0.0 };
  for (std::size_t b = 0; b < 4; ++b)
  {
    dense[b][1] = static_cast<double>(100 + b);
  }
  const VectorMatrix<double, 4>& cdense = dense;
  VectorMatrix<double, 4>::ConstGroupView dcgv(cdense, 0);
  auto dv = dcgv.GetConstColumnView(1);
  static_assert(GroupedDenseMatrixColumnView<decltype(dv)>);
  for (std::size_t i = 0; i < 4; ++i)
  {
    EXPECT_EQ(dv.base_[i], static_cast<double>(100 + i));
  }
}

// ---------------------------------------------------------------------------
// GroupView.Fill / GroupView.Copy: bulk-fill and bulk-copy sibling primitives
// of ForEachBlock. Each test verifies:
//   (a) the write only touches the addressed block (surrounding non-zero
//       elements and neighbouring blocks in the group are untouched);
//   (b) semantically the result matches the ForEachBlock equivalent.
// ---------------------------------------------------------------------------
namespace
{
  template<class SM, std::size_t L>
  void RunSparseFillCopyTest()
  {
    // Two non-zero positions per block; a full group and one partial group.
    auto builder = SM::Create(3).WithElement(0, 1).WithElement(2, 0).SetNumberOfBlocks(L + 1);
    SM matrix{ builder };
    for (std::size_t b = 0; b < L + 1; ++b)
    {
      matrix[b][0][1] = static_cast<double>(b + 1);
      matrix[b][2][0] = static_cast<double>(100 + b);
    }

    // Fill on the full group: only [0][1] cells change, [2][0] cells untouched.
    {
      typename SM::GroupView gv(matrix, 0);
      gv.Fill(gv.GetBlockView(matrix.VectorIndex(0, 0, 1)), 42.0);
      for (std::size_t b = 0; b < L; ++b)
      {
        EXPECT_DOUBLE_EQ(matrix[b][0][1], 42.0);
        EXPECT_DOUBLE_EQ(matrix[b][2][0], static_cast<double>(100 + b));
      }
      // Partial group is unaffected.
      EXPECT_DOUBLE_EQ(matrix[L][0][1], static_cast<double>(L + 1));
    }

    // Copy on the full group: [2][0] <- [0][1].
    {
      typename SM::GroupView gv(matrix, 0);
      gv.Copy(gv.GetBlockView(matrix.VectorIndex(0, 2, 0)), gv.GetConstBlockView(matrix.VectorIndex(0, 0, 1)));
      for (std::size_t b = 0; b < L; ++b)
      {
        EXPECT_DOUBLE_EQ(matrix[b][2][0], 42.0);
        EXPECT_DOUBLE_EQ(matrix[b][0][1], 42.0);
      }
    }

    // Fill on the partial (last) group behaves the same as ForEachBlock
    // semantically for the in-range cells. Because the storage is padded to L,
    // padding cells past num_blocks_in_group are also written but never
    // observable through the matrix interface.
    {
      const std::size_t last_group = matrix.NumberOfBlocks() / L;
      typename SM::GroupView gv(matrix, last_group);
      gv.Fill(gv.GetBlockView(matrix.VectorIndex(0, 0, 1)), -7.0);
      EXPECT_DOUBLE_EQ(matrix[L][0][1], -7.0);
    }
  }
}  // namespace

TEST(GroupedView, FillAndCopySparseStandardCSR)
{
  using SM = SparseMatrix<double, SparseMatrixStandardOrderingCompressedSparseRow>;
  // Standard ordering: L=1 (one "group" per block).
  auto builder = SM::Create(3).WithElement(0, 1).WithElement(2, 0).SetNumberOfBlocks(2);
  SM matrix{ builder };
  for (std::size_t b = 0; b < 2; ++b)
  {
    matrix[b][0][1] = static_cast<double>(b + 1);
    matrix[b][2][0] = static_cast<double>(100 + b);
  }

  SM::GroupView gv0(matrix, 0);
  gv0.Fill(gv0.GetBlockView(matrix.VectorIndex(0, 0, 1)), 42.0);
  EXPECT_DOUBLE_EQ(matrix[0][0][1], 42.0);
  EXPECT_DOUBLE_EQ(matrix[0][2][0], 100.0);  // untouched
  EXPECT_DOUBLE_EQ(matrix[1][0][1], 2.0);    // other group untouched

  gv0.Copy(gv0.GetBlockView(matrix.VectorIndex(0, 2, 0)), gv0.GetConstBlockView(matrix.VectorIndex(0, 0, 1)));
  EXPECT_DOUBLE_EQ(matrix[0][2][0], 42.0);
}

TEST(GroupedView, FillAndCopySparseStandardCSC)
{
  using SM = SparseMatrix<double, SparseMatrixStandardOrderingCompressedSparseColumn>;
  auto builder = SM::Create(3).WithElement(0, 1).WithElement(2, 0).SetNumberOfBlocks(2);
  SM matrix{ builder };
  for (std::size_t b = 0; b < 2; ++b)
  {
    matrix[b][0][1] = static_cast<double>(b + 1);
    matrix[b][2][0] = static_cast<double>(100 + b);
  }

  SM::GroupView gv1(matrix, 1);
  gv1.Fill(gv1.GetBlockView(matrix.VectorIndex(0, 0, 1)), 9.0);
  EXPECT_DOUBLE_EQ(matrix[1][0][1], 9.0);
  EXPECT_DOUBLE_EQ(matrix[0][0][1], 1.0);  // other group untouched
  EXPECT_DOUBLE_EQ(matrix[1][2][0], 101.0);

  gv1.Copy(gv1.GetBlockView(matrix.VectorIndex(0, 2, 0)), gv1.GetConstBlockView(matrix.VectorIndex(0, 0, 1)));
  EXPECT_DOUBLE_EQ(matrix[1][2][0], 9.0);
}

TEST(GroupedView, FillAndCopySparseVectorCSR_L4)
{
  RunSparseFillCopyTest<SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<4>>, 4>();
}

TEST(GroupedView, FillAndCopySparseVectorCSC_L4)
{
  RunSparseFillCopyTest<SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseColumn<4>>, 4>();
}

TEST(GroupedView, FillAndCopySparseVectorCSR_L8)
{
  RunSparseFillCopyTest<SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<8>>, 8>();
}

// Fill/Copy on a mutable GroupView must accept a GroupedConstBlockView from a
// separate ConstGroupView on the same underlying storage as the Copy source.
TEST(GroupedView, FillCopyCrossGroupViewSources)
{
  using SM = SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  auto builder = SM::Create(2).WithElement(0, 1).SetNumberOfBlocks(4);
  SM src{ builder };
  SM dst{ builder };
  for (std::size_t b = 0; b < 4; ++b)
  {
    src[b][0][1] = static_cast<double>(b + 1);
    dst[b][0][1] = 0.0;
  }

  SM::GroupView dgv(dst, 0);
  SM::ConstGroupView sgv(src, 0);
  dgv.Copy(dgv.GetBlockView(dst.VectorIndex(0, 0, 1)), sgv.GetConstBlockView(src.VectorIndex(0, 0, 1)));
  for (std::size_t b = 0; b < 4; ++b)
  {
    EXPECT_DOUBLE_EQ(dst[b][0][1], static_cast<double>(b + 1));
  }

  // Fill result matches the equivalent ForEachBlock lambda.
  SM ref{ builder };
  for (std::size_t b = 0; b < 4; ++b)
  {
    ref[b][0][1] = 999.0;
  }
  SM::GroupView rgv(ref, 0);
  rgv.ForEachBlock([](double& x) { x = 55.5; }, rgv.GetBlockView(ref.VectorIndex(0, 0, 1)));
  dgv.Fill(dgv.GetBlockView(dst.VectorIndex(0, 0, 1)), 55.5);
  for (std::size_t b = 0; b < 4; ++b)
  {
    EXPECT_DOUBLE_EQ(dst[b][0][1], ref[b][0][1]);
  }
}
