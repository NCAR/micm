// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/phase.hpp>

#include <cstddef>
#include <string>

namespace micm
{
  struct SurfaceRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Gas-phase species reacting on surface
    PhaseSpecies phase_species_;
    /// @brief Reaction probability (0-1) [unitless]
    double reaction_probability_{ 1.0 };
  };

  /// @brief GPU-safe calculation data for a surface reaction.
  ///        Populated by ReactionRateConstantStore::BuildFrom from SurfaceRateConstantParameters;
  ///        do not construct directly.
  struct SurfaceRateConstantData
  {
    /// @brief Gas-phase diffusion coefficient for the reacting species [m2 s-1]
    double diffusion_coefficient_;
    /// @brief Precomputed factor for mean free speed: 8 * R / (pi * Mw) [K-1 m2 s-2]
    double mean_free_speed_factor_;
    /// @brief Reaction probability (0-1) [unitless]
    double reaction_probability_;
    /// @brief Index into custom_rate_parameters_[cell] for aerosol effective radius [m];
    ///        particle number concentration [# m-3] is at custom_param_base_index_ + 1
    std::size_t custom_param_base_index_;
  };
}  // namespace micm
