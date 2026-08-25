// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <Kokkos_Core.hpp>
#include <algorithm>

namespace micm::detail
{
  /// @brief Team size for a launch whose intra-team loop is TeamThreadRange(team, L)
  /// @tparam L Trip count of the intra-team loop
  /// @param team_functor The functor the team policy will be launched with
  /// @return A team size to hand to Kokkos::TeamPolicy
  template<Index L, typename Functor>
  int TeamSizeForL(const Functor& team_functor)
  {
    static const int team_size = [&team_functor]()
    {
      const int max_team_size = Kokkos::TeamPolicy<>(1, Kokkos::AUTO).team_size_max(team_functor, Kokkos::ParallelForTag());
      // A team size of 0 is not a legal TeamPolicy argument
      return std::max(1, std::min(static_cast<int>(L), max_team_size));
    }();
    return team_size;
  }
}  // namespace micm::detail
