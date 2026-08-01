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
    Kokkos::Array<T, L> storage_;

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

}  // namespace micm
