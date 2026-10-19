// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <type_traits>

namespace micm::detail
{
  /// @brief Team size for a launch whose intra-team loop is TeamThreadRange(team, L)
  /// @tparam L Trip count of the intra-team loop
  /// @param team_functor The functor the team policy will be launched with
  /// @return A team size to hand to Kokkos::TeamPolicy
  template<Index L, typename Functor>
  int TeamSizeForL(const Functor& team_functor)
  {
    // On a host execution space (OpenMP, Threads, Serial), each team member is a CPU thread, and the
    // intra-team loops end with a team barrier. A team of several threads on one group of L cells meets at
    // a barrier after every short row operation, and the host barrier sleeps when it waits, so more threads
    // made the solve much slower. One thread for each team runs whole groups on each thread instead,
    // and the compiler vectorizes the loop over L.
    if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>)
    {
      return 1;
    }
    static const int team_size = [&team_functor]()
    {
      const int max_team_size = Kokkos::TeamPolicy<>(1, Kokkos::AUTO).team_size_max(team_functor, Kokkos::ParallelForTag());
      // A team size of 0 is not a legal TeamPolicy argument
      return std::max(1, std::min(static_cast<int>(L), max_team_size));
    }();
    return team_size;
  }
}  // namespace micm::detail
