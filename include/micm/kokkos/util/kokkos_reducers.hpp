// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_scalar_view.hpp>

#include <Kokkos_Core.hpp>

namespace micm
{
  /// @brief Sum reduction (`acc += x`).
  template<typename T>
  struct KokkosSum
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using value_type = T;
    Kokkos::View<T*> view_;

    KOKKOS_INLINE_FUNCTION constexpr explicit KokkosSum(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const  // NOLINT(readability-identifier-naming)
    {
      return view_.data();
    }

    KOKKOS_INLINE_FUNCTION static constexpr T Identity()
    {
      return T{};
    }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const
    {
      T local = Identity();
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::Sum<T>(local));
      Kokkos::single(Kokkos::PerTeam(team), [&]() { Kokkos::atomic_add(view_.data(), local); });
    }

    KOKKOS_INLINE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      dst += src;
    }
  };

  /// @brief Max reduction (`acc = std::max(acc, x)`).
  template<typename T>
  struct KokkosMax
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using value_type = T;
    Kokkos::View<T*> view_;

    KOKKOS_INLINE_FUNCTION constexpr explicit KokkosMax(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const  // NOLINT(readability-identifier-naming)
    {
      return view_.data();
    }

    KOKKOS_INLINE_FUNCTION static constexpr T Identity()
    {
      return T{};
    }
    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const
    {
      T local = Identity();
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::Max<T>(local));
      Kokkos::single(Kokkos::PerTeam(team), [&]() { Kokkos::atomic_fetch_max(view_.data(), local); });
    }

    KOKKOS_INLINE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      if (src > dst)
      {
        dst = src;
      }
    }
  };

  /// @brief Logical And redicer (`acc = acc || x`).
  struct KokkosLOr
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using T = Bool;
    using value_type = T;
    Kokkos::View<T*> view_;

    KOKKOS_INLINE_FUNCTION constexpr explicit KokkosLOr(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const
    {
      return view_.data();
    }

    KOKKOS_INLINE_FUNCTION static constexpr T Identity()
    {
      return T{};
    }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const
    {
      Bool local = Identity();
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::LOr<Bool>(local));
      Kokkos::single(
          Kokkos::PerTeam(team),
          [&]()
          {
            if (local)
              Kokkos::atomic_store(view_.data(), true);
          });
    }

    KOKKOS_INLINE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      dst = dst || src;
    }
  };

  /// @brief Logical Or reducer (`acc = acc && x`).
  struct KokkosLAnd
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using T = Bool;
    using value_type = T;
    Kokkos::View<T*> view_;

    KOKKOS_INLINE_FUNCTION constexpr explicit KokkosLAnd(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const
    {
      return view_.data();
    }

    KOKKOS_INLINE_FUNCTION static constexpr T Identity()
    {
      return T{};
    }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const
    {
      Bool local = Identity();
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::LAnd<Bool>(local));
      Kokkos::single(
          Kokkos::PerTeam(team),
          [&]()
          {
            if (!local)
              Kokkos::atomic_store(view_.data(), false);
          });
    }

    KOKKOS_INLINE_FUNCTION static constexpr void Join(T& dst, const T& src)
    {
      dst = dst && src;
    }
  };
}  // namespace micm