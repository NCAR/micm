// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/rate_constant/rate_constant.hpp>
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

    Yield(const Species& species)
        : species_(species),
          coefficient_(1.0)
    {
    }
  };

}  // namespace micm