// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

namespace micm
{
  struct BranchedRateConstantParameters
  {
    enum class Branch
    {
      Alkoxy,
      Nitrate
    };
    /// @brief reaction branch
    Branch branch_;
    /// @brief pre-exponential factor
    Real X_;
    /// @brief exponential factor
    Real Y_;
    /// @brief branching factor
    Real a0_;
    /// @brief number of heavy atoms in the RO2 reacting species (excluding the peroxy moiety)
    Index n_;
    /// @brief Precomputed low-pressure rate factor: 2e-22 * N_A * 1e-6 * exp(n_)
    ///        Set by ReactionRateConstantStore::BuildFrom; do not set manually.
    Real k0_{ 0.0 };
    /// @brief Precomputed branching ratio factor: A(293, [M]_ref) * (1 - a0_) / a0_
    ///        Set by ReactionRateConstantStore::BuildFrom; do not set manually.
    Real z_{ 0.0 };
  };
}  // namespace micm
