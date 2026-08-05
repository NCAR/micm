// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_scalar_view.hpp>

namespace micm
{
  /// @brief Sum reduction (`acc += x`).
  template<typename T>
  struct KokkosSum
  {
    using value_type = T;
    Kokkos::View<T*> view_;

    constexpr explicit KokkosSum(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const { return view_.data(); }

    KOKKOS_INLINE_FUNCTION static constexpr T identity() { return T{}; }
  };
}