// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <Kokkos_Core.hpp>
#include <mutex>

namespace micm
{
  namespace kokkos
  {
    /// @brief Initialize Kokkos if not already initialized
    inline void Initialize()
    {
      static std::once_flag init_flag;
      std::call_once(init_flag, []() {
        if (!Kokkos::is_initialized())
        {
          Kokkos::initialize();
        }
      });
    }

    /// @brief Finalize Kokkos if initialized
    inline void Finalize()
    {
      if (Kokkos::is_initialized())
      {
        Kokkos::finalize();
      }
    }

    /// @brief Scope guard for Kokkos initialization and finalization
    struct ScopeGuard
    {
      ScopeGuard(int& argc, char* argv[])
      {
        Kokkos::initialize(argc, argv);
      }
      ~ScopeGuard()
      {
        Kokkos::finalize();
      }
    };
  }  // namespace kokkos
}  // namespace micm
