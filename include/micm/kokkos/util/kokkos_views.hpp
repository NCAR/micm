// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <Kokkos_Core.hpp>

namespace micm
{

  /// @brief Device-callable, ungrouped column view for a Kokkos-backed dense matrix.
  ///
  /// Mirrors `VectorMatrix::ColumnView`, but stores a raw pointer directly into the
  /// matrix's flat storage (the `Kokkos::View`'s data pointer) instead of a pointer
  /// back to the host matrix object. This makes it trivially copyable into a
  /// `KOKKOS_LAMBDA`, where dereferencing a host object pointer would be unsafe.
  ///
  /// `L` is the row-group size (`VectorMatrix`'s tiered grouping factor). Given a
  /// group index `g` and `row_in_group`, the corresponding element lives at
  /// `data_[(g * y_dim_ + column_index_) * L + row_in_group]` -- the same layout
  /// `VectorMatrix` uses.
  template<class T, Index L>
  class KokkosColumnView
  {
    T* data_;
    Index y_dim_;
    Index column_index_;

   public:
    using category = DenseMatrixColumnViewTag;

    KOKKOS_INLINE_FUNCTION
    KokkosColumnView(T* data, Index y_dim, Index column_index)
        : data_(data),
          y_dim_(y_dim),
          column_index_(column_index)
    {
    }

    KOKKOS_INLINE_FUNCTION Index ColumnIndex() const
    {
      return column_index_;
    }

    KOKKOS_INLINE_FUNCTION T* Data() const
    {
      return data_;
    }

    KOKKOS_INLINE_FUNCTION Index YDim() const
    {
      return y_dim_;
    }
  };

  /// @brief Const variant of KokkosColumnView. See KokkosColumnView for details.
  template<class T, Index L>
  class KokkosConstColumnView
  {
    const T* data_;
    Index y_dim_;
    Index column_index_;

   public:
    using category = DenseMatrixColumnViewTag;

    KOKKOS_INLINE_FUNCTION
    KokkosConstColumnView(const T* data, Index y_dim, Index column_index)
        : data_(data),
          y_dim_(y_dim),
          column_index_(column_index)
    {
    }

    KOKKOS_INLINE_FUNCTION Index ColumnIndex() const
    {
      return column_index_;
    }

    KOKKOS_INLINE_FUNCTION const T* Data() const
    {
      return data_;
    }

    KOKKOS_INLINE_FUNCTION Index YDim() const
    {
      return y_dim_;
    }
  };

  /// @brief Enriched mutable column view for a single row-group of a Kokkos-backed
  ///        dense matrix.
  ///
  /// Carries a precomputed `base_` pointer at the first row of the group's L-row
  /// block for a given column, so element access within the group is the
  /// contiguous `base_[row_in_group]` instead of recomputing
  /// `(group * y_dim + column) * L + row_in_group` per element. Mirrors
  /// `VectorMatrix::GroupView::GroupedColumnView`.
  template<class T>
  struct KokkosGroupedColumnView
  {
    using category = GroupedDenseMatrixColumnViewTag;
    T* base_;
  };

  /// @brief Const variant of KokkosGroupedColumnView. See KokkosGroupedColumnView
  ///        for details.
  template<class T>
  struct KokkosGroupedConstColumnView
  {
    using category = GroupedDenseMatrixColumnViewTag;
    const T* base_;
  };

  /// @brief A row-local temporary variable with its own device-callable storage.
  ///
  /// Mirrors `VectorMatrix::RowVariable`, but uses `Kokkos::Array` rather than
  /// `std::array` for the backing storage, and marks its accessors
  /// `KOKKOS_INLINE_FUNCTION` so instances can be constructed and used inside a
  /// `KOKKOS_LAMBDA` (including on GPU backends).
  template<class T, Index L>
  class KokkosRowVariable
  {
    Kokkos::Array<T, L> storage_;  // TODO look at KokkosBlockVariable below and see if we can use that here

   public:
    using category = BlockVariableTag;

    KOKKOS_DEFAULTED_FUNCTION
    KokkosRowVariable() = default;

    KOKKOS_INLINE_FUNCTION Kokkos::Array<T, L>& Get()
    {
      return storage_;
    }

    KOKKOS_INLINE_FUNCTION const Kokkos::Array<T, L>& Get() const
    {
      return storage_;
    }
  };

  /// @brief Device-callable, ungrouped block view for a Kokkos-backed sparse matrix.
  ///
  /// Mirrors `SparseMatrix::BlockView`/`ConstBlockView`, but stores a raw pointer
  /// directly into the matrix's flat storage (the `Kokkos::View`'s data pointer)
  /// instead of a pointer back to the host matrix object, so it can be captured by
  /// value into a `KOKKOS_LAMBDA`.
  ///
  /// `L` is the block-group size (the sparse ordering policy's `GroupVectorSize()`).
  /// Given a block index `b`, the corresponding element lives at
  /// `data_[(b / L) * flat_block_size_ * L + element_position_ + b % L]`, where
  /// `element_position_` is the block-relative offset for this (row, column)
  /// returned by the ordering policy (the same value `SparseMatrix::ElementPosition()`
  /// exposes on the host).
  template<class T, Index L>
  class KokkosBlockView
  {
    T* data_;
    Index flat_block_size_;
    Index element_position_;

   public:
    using category = SparseMatrixBlockViewTag;

    KOKKOS_INLINE_FUNCTION
    KokkosBlockView(T* data, Index flat_block_size, Index element_position)
        : data_(data),
          flat_block_size_(flat_block_size),
          element_position_(element_position)
    {
    }

    KOKKOS_INLINE_FUNCTION Index ElementPosition() const
    {
      return element_position_;
    }

    KOKKOS_INLINE_FUNCTION T* Data() const
    {
      return data_;
    }

    KOKKOS_INLINE_FUNCTION Index FlatBlockSize() const
    {
      return flat_block_size_;
    }
  };

  /// @brief Const variant of KokkosBlockView. See KokkosBlockView for details.
  template<class T, Index L>
  class KokkosConstBlockView
  {
    const T* data_;
    Index flat_block_size_;
    Index element_position_;

   public:
    using category = SparseMatrixBlockViewTag;

    KOKKOS_INLINE_FUNCTION
    KokkosConstBlockView(const T* data, Index flat_block_size, Index element_position)
        : data_(data),
          flat_block_size_(flat_block_size),
          element_position_(element_position)
    {
    }

    KOKKOS_INLINE_FUNCTION Index ElementPosition() const
    {
      return element_position_;
    }

    KOKKOS_INLINE_FUNCTION const T* Data() const
    {
      return data_;
    }

    KOKKOS_INLINE_FUNCTION Index FlatBlockSize() const
    {
      return flat_block_size_;
    }
  };

  /// @brief Enriched mutable block view for a single block-group of a Kokkos-backed
  ///        sparse matrix.
  ///
  /// Carries a precomputed `group_base_` pointer at the start of the group's slice
  /// of the sparse data vector, so element access within the group is the
  /// contiguous `group_base_[block_offset_ + block_in_group]` instead of
  /// recomputing the flat index per element. Mirrors the ordering policy's
  /// `GroupView::GroupedBlockView`.
  template<class T>
  struct KokkosGroupedBlockView
  {
    using category = GroupedSparseMatrixBlockViewTag;
    T* group_base_;
    Index block_offset_;
  };

  /// @brief Const variant of KokkosGroupedBlockView. See KokkosGroupedBlockView for
  ///        details.
  template<class T>
  struct KokkosGroupedConstBlockView
  {
    using category = GroupedSparseMatrixBlockViewTag;
    const T* group_base_;
    Index block_offset_;
  };

  /// @brief A block-local temporary variable with its own device-callable storage.
  ///
  /// Mirrors the ordering policy's `BlockVariable`: a `Kokkos::Array<T, L>` when
  /// blocks are grouped (`L > 1`), or a bare scalar `T` for standard ordering
  /// (`L == 1`). Accessors are `KOKKOS_INLINE_FUNCTION` so instances can be
  /// constructed and used inside a `KOKKOS_LAMBDA`.
  template<class T, Index L>
  class KokkosBlockVariable
  {
    std::conditional_t<(L > 1), Kokkos::Array<T, L>, T> storage_;

   public:
    using category = BlockVariableTag;

    KOKKOS_DEFAULTED_FUNCTION
    KokkosBlockVariable() = default;

    KOKKOS_INLINE_FUNCTION auto& Get()
    {
      return storage_;
    }

    KOKKOS_INLINE_FUNCTION const auto& Get() const
    {
      return storage_;
    }
  };

}  // namespace micm
