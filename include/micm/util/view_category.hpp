// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <concepts>
#include <type_traits>

namespace micm
{

  // ============================================================================
  // View Category Tags - What kind of view
  // ============================================================================

  /// @brief Tag for sparse matrix block views (have RowIndex + ColumnIndex)
  struct SparseMatrixBlockViewTag
  {
  };

  /// @brief Tag for dense matrix column views (have ColumnIndex + GetMatrix)
  struct DenseMatrixColumnViewTag
  {
  };

  /// @brief Tag for block variables (vector-like data holders)
  struct BlockVariableTag
  {
  };

  // ============================================================================
  // Grouping Strategy Tags - How data is organized
  // ============================================================================

  /// @brief Simple grouping: L==1, group index directly maps to element
  /// Used by: Matrix (always), VectorMatrix (when L==1),
  ///          Standard ordering sparse (always), Vector ordering sparse (when L==1)
  struct SimpleGroupingTag
  {
  };

  /// @brief Tiered grouping: L>1, groups contain L elements, need block_in_group offset
  /// Used by: VectorMatrix (when L>1), Vector ordering sparse (when L>1)
  struct TieredGroupingTag
  {
  };

  // ============================================================================
  // View Category Trait
  // ============================================================================

  /// @brief Helper to check if a type has a nested 'category' type
  template<typename T, typename = void>
  struct has_category : std::false_type
  {
  };

  template<typename T>
  struct has_category<T, std::void_t<typename T::category>> : std::true_type
  {
  };

  /// @brief Determines the category of a view type (checks for nested 'category' type first)
  template<typename T, typename = void>
  struct ViewCategory;

  /// @brief If type has a nested 'category' type, use it
  template<typename T>
  struct ViewCategory<T, std::enable_if_t<has_category<std::remove_cvref_t<T>>::value>>
  {
    using type = typename std::remove_cvref_t<T>::category;
  };

  /// @brief Helper alias
  template<typename T>
  using ViewCategory_t = typename ViewCategory<std::remove_cvref_t<T>>::type;

  // ============================================================================
  // Grouping Strategy Trait
  // ============================================================================

  /// @brief Determines the grouping strategy of a matrix type (no default - must be specialized)
  template<typename T>
  struct GroupingStrategy;

  /// @brief Helper alias
  template<typename T>
  using GroupingStrategy_t = typename GroupingStrategy<std::remove_cvref_t<T>>::type;

  // ============================================================================
  // Concepts for View Categories
  // ============================================================================

  /// @brief Sparse matrix block view (has RowIndex/ColumnIndex)
  template<typename T>
  concept SparseMatrixBlockView = std::same_as<ViewCategory_t<T>, SparseMatrixBlockViewTag>;

  /// @brief Dense matrix column view (has ColumnIndex/GetMatrix)
  template<typename T>
  concept DenseMatrixColumnView = std::same_as<ViewCategory_t<T>, DenseMatrixColumnViewTag>;

  /// @brief Block variable (vector-like data holder)
  template<typename T>
  concept BlockVariableView = std::same_as<ViewCategory_t<T>, BlockVariableTag>;

  /// @brief Vector-like type (has operator[] and size())
  /// Excludes Matrix types, view types, and proxy types to avoid ambiguity
  template<typename T>
  concept VectorLike = requires(T t, std::size_t i) {
    { t[i] };  // Can index with []
    { t.size() } -> std::convertible_to<std::size_t>;
  } && !DenseMatrixColumnView<T> && !BlockVariableView<T> && !SparseMatrixBlockView<T> &&
       !requires(T t) { t.NumRows(); t.NumColumns(); };  // Exclude matrix types

  // ============================================================================
  // Concepts for Grouping Strategies
  // ============================================================================

  /// @brief Matrix uses simple grouping (L==1)
  template<typename T>
  concept HasSimpleGrouping = requires { typename GroupingStrategy<std::remove_cvref_t<T>>::type; } &&
                              std::same_as<GroupingStrategy_t<T>, SimpleGroupingTag>;

  /// @brief Matrix uses tiered grouping (L>1)
  template<typename T>
  concept HasTieredGrouping = requires { typename GroupingStrategy<std::remove_cvref_t<T>>::type; } &&
                              std::same_as<GroupingStrategy_t<T>, TieredGroupingTag>;

}  // namespace micm
