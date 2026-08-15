// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <algorithm>
#include <limits>

namespace micm
{
  /// @file
  /// @brief Reducer types for `Reduce`/`ReduceStrict` on matrix policy GroupViews.
  ///
  /// These mirror the shape of `Kokkos::Sum` / `Kokkos::Max` / `Kokkos::LOr` /
  /// `Kokkos::LAnd`: each reducer wraps a caller-owned destination scalar (`T&`)
  /// and exposes:
  /// - `value_type`: the accumulator's scalar type.
  /// - `Reference()`: returns the destination scalar by reference.
  /// - `Identity()`: the identity value used to seed a fresh accumulator.
  /// - `Join(dst, src)`: how per-thread accumulators combine into `dst`.
  ///
  /// On host matrix policies, `Reduce` iterates serially and writes into
  /// `Reference()` directly. On `KokkosDenseMatrix`, the reducer is translated
  /// into the matching `Kokkos::Sum` / `Kokkos::Max` / `Kokkos::LOr` /
  /// `Kokkos::LAnd` before dispatching `Kokkos::parallel_reduce`.

  /// @brief Sum reduction (`acc += x`).
  template<typename T>
  struct Sum
  {
    using value_type = T;
    T& value_;

    MICM_DEVICE_FUNCTION constexpr explicit Sum(T& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    MICM_DEVICE_FUNCTION explicit Sum(const Scalar& scalar)
        : value_(scalar.HostValue())
    {
    }

    MICM_DEVICE_FUNCTION constexpr T& Reference() const
    {
      return value_;
    }

    MICM_DEVICE_FUNCTION static constexpr T Identity()
    {
      return T{};
    }

    MICM_DEVICE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      dst += src;
    }
  };

  template<typename T>
  Sum(T&) -> Sum<T>;

  template<typename Scalar>
  Sum(const Scalar&) -> Sum<typename Scalar::value_type>;

  /// @brief Max reduction (`acc = max(acc, x)`).
  template<typename T>
  struct Max
  {
    using value_type = T;
    T& value_;

    MICM_DEVICE_FUNCTION constexpr explicit Max(T& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    MICM_DEVICE_FUNCTION explicit Max(const Scalar& scalar)
        : value_(scalar.HostValue())
    {
    }

    MICM_DEVICE_FUNCTION constexpr T& Reference() const
    {
      return value_;
    }

    MICM_DEVICE_FUNCTION static constexpr T Identity()
    {
      return std::numeric_limits<T>::lowest();
    }

    MICM_DEVICE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      if (src > dst)
      {
        dst = src;
      }
    }
  };

  template<typename T>
  Max(T&) -> Max<T>;

  template<typename Scalar>
  Max(const Scalar&) -> Max<typename Scalar::value_type>;

  /// @brief Logical-OR reduction (`acc = acc || x`).
  struct LOr
  {
    using value_type = Bool;
    Bool& value_;

    MICM_DEVICE_FUNCTION constexpr explicit LOr(Bool& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    MICM_DEVICE_FUNCTION explicit LOr(const Scalar& scalar)
        : value_(scalar.HostValue())
    {
    }

    MICM_DEVICE_FUNCTION constexpr Bool& Reference() const
    {
      return value_;
    }

    MICM_DEVICE_FUNCTION static constexpr Bool Identity()
    {
      return static_cast<Bool>(false);
    }

    MICM_DEVICE_FUNCTION static constexpr void Join(Bool& dst, Bool src)
    {
      dst = static_cast<Bool>(static_cast<bool>(dst) || static_cast<bool>(src));
    }
  };

  /// @brief Logical-AND reduction (`acc = acc && x`).
  struct LAnd
  {
    using value_type = Bool;
    Bool& value_;

    MICM_DEVICE_FUNCTION constexpr explicit LAnd(Bool& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    MICM_DEVICE_FUNCTION explicit LAnd(const Scalar& scalar)
        : value_(scalar.HostValue())
    {
    }

    MICM_DEVICE_FUNCTION constexpr Bool& Reference() const
    {
      return value_;
    }

    MICM_DEVICE_FUNCTION static constexpr Bool Identity()
    {
      return static_cast<Bool>(true);
    }

    MICM_DEVICE_FUNCTION static constexpr void Join(Bool& dst, Bool src)
    {
      dst = static_cast<Bool>(static_cast<bool>(dst) && static_cast<bool>(src));
    }
  };
}  // namespace micm
