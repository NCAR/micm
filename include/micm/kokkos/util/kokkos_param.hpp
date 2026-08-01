// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <Kokkos_Core.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace micm
{
  namespace kokkos
  {
    /// @brief Copy a host std::vector into a new Kokkos::View on the default execution space
    template<typename T>
    Kokkos::View<T*> CopyVectorToView(const std::string& name, const std::vector<T>& src)
    {
      Kokkos::View<T*> view(name, src.size());
      auto host = Kokkos::create_mirror_view(view);
      for (std::size_t i = 0; i < src.size(); ++i)
      {
        host(i) = src[i];
      }
      Kokkos::deep_copy(view, host);
      return view;
    }
  }  // namespace kokkos
}  // namespace micm
