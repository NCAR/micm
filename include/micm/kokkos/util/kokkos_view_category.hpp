// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/view_category.hpp>

namespace micm
{
  // ============================================================================
  // Concepts for View Categories
  // ============================================================================

  /// @brief Vector-like with padded cells and device data accessibility
  template<typename T>
  concept KokkosVectorLike = VectorLike<T> && requires(T t) {
    { t.GetView() };
  };
}  // namespace micm