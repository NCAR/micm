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

    constexpr explicit KokkosSum(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const { return view_.data(); }

    KOKKOS_INLINE_FUNCTION static constexpr T identity() { return T{}; }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const {
        T local = identity();
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::Sum<T>(local));
        Kokkos::single(Kokkos::PerTeam(team), [&]() { Kokkos::atomic_add(view_.data(), local); });
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

    constexpr explicit KokkosMax(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const { return view_.data(); }

    KOKKOS_INLINE_FUNCTION static constexpr T identity() { return T{}; }
    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const {
        T local = identity();
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::Max<T>(local));
        Kokkos::single(Kokkos::PerTeam(team), [&]() { Kokkos::atomic_fetch_max(view_.data(), local); });
    }
  };

  /// @brief Logical And redicer (`acc = acc || x`).
  struct KokkosLOr
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using T = bool;
    using value_type = T;
    Kokkos::View<T*> view_;

    constexpr explicit KokkosLOr(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const { return view_.data(); }

    KOKKOS_INLINE_FUNCTION static constexpr T identity() { return T{}; }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const {
        bool local = identity();
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::LOr<bool>(local));
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            if (local) Kokkos::atomic_store(view_.data(), true);
        });
    }
  };

  /// @brief Logical Or reducer (`acc = acc && x`).
  struct KokkosLAnd
  {
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    using T = bool;
    using value_type = T;
    Kokkos::View<T*> view_;

    constexpr explicit KokkosLAnd(const KokkosScalarView<T>& scalar)
        : view_(scalar.GetDeviceView())
    {
    }

    KOKKOS_INLINE_FUNCTION T* device_ptr() const { return view_.data(); }

    KOKKOS_INLINE_FUNCTION static constexpr T identity() { return T{}; }

    template<typename RowFunc>
    KOKKOS_INLINE_FUNCTION void team_reduce(const TeamMember& team, Index count, RowFunc&& row_func) const {
        bool local = identity();
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, count), row_func, Kokkos::LAnd<bool>(local));
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            if (!local) Kokkos::atomic_store(view_.data(), false);
        });
    }
  };
}