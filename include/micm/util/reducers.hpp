// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

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
  /// - `reference()`: returns the destination scalar by reference.
  /// - `identity()`: the identity value used to seed a fresh accumulator.
  /// - `join(dst, src)`: how per-thread accumulators combine into `dst`.
  ///
  /// On host matrix policies, `Reduce` iterates serially and writes into
  /// `reference()` directly. On `KokkosDenseMatrix`, the reducer is translated
  /// into the matching `Kokkos::Sum` / `Kokkos::Max` / `Kokkos::LOr` /
  /// `Kokkos::LAnd` before dispatching `Kokkos::parallel_reduce`.

  /// @brief Sum reduction (`acc += x`).
  template<typename T>
  struct Sum
  {
    using value_type = T;
    T& value_;

    constexpr explicit Sum(T& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    explicit Sum(const Scalar& scalar)
        : value_(scalar.host_value())
    {
    }

    constexpr T& reference() const
    {
      return value_;
    }

    static constexpr T identity()
    {
      return T{};
    }

    static constexpr void join(T& dst, const T& src)
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

    constexpr explicit Max(T& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    explicit Max(const Scalar& scalar)
        : value_(scalar.host_value())
    {
    }

    constexpr T& reference() const
    {
      return value_;
    }

    static constexpr T identity()
    {
      return std::numeric_limits<T>::lowest();
    }

    static constexpr void join(T& dst, const T& src)
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
    using value_type = bool;
    bool& value_;

    constexpr explicit LOr(bool& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    explicit LOr(const Scalar& scalar)
        : value_(scalar.host_value())
    {
    }

    constexpr bool& reference() const
    {
      return value_;
    }

    static constexpr bool identity()
    {
      return false;
    }

    static constexpr void join(bool& dst, bool src)
    {
      dst = dst || src;
    }
  };

  /// @brief Logical-AND reduction (`acc = acc && x`).
  struct LAnd
  {
    using value_type = bool;
    bool& value_;

    constexpr explicit LAnd(bool& value)
        : value_(value)
    {
    }

    template<typename Scalar>
    explicit LAnd(const Scalar& scalar)
        : value_(scalar.host_value())
    {
    }

    constexpr bool& reference() const
    {
      return value_;
    }

    static constexpr bool identity()
    {
      return true;
    }

    static constexpr void join(bool& dst, bool src)
    {
      dst = dst && src;
    }
  };
}  // namespace micm
