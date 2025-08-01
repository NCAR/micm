// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <memory>
#include <vector>

namespace micm
{
  /// @brief Represents a product in a chemical reaction, associating a species with
  ///        its stoichiometric coefficient
  struct Yield
  {
    Species species_;
    double coefficient_;

    Yield() = default;

    Yield(const Species& species, double coefficient)
        : species_(species),
          coefficient_(coefficient)
    {
    }
  };

  /// @brief Represents a homogeneous or single-phase reaction with phase-specific
  ///        reactants, products, and a rate constant
  struct PhaseReaction
  {
    Phase phase_;
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;
  };

  /// @brief Represents a phase-transfer reaction where reactants move from one phase to another,
  ///        producing products in the destination phase
  struct PhaseTransferReaction
  {
    Phase source_phase_;
    Phase destination_phase_;
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;
  };

}  // namespace micm