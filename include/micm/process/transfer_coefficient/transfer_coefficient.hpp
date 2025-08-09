// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/conditions.hpp>

#include <memory>
#include <utility>


namespace micm
{

  class TransferCoefficient
  {
  public:
    virtual ~TransferCoefficient() = default;

    virtual std::unique_ptr<TransferCoefficient> Clone() const = 0;

    virtual double Calculate() const
    {
      return 0.0;
    }

    /// @brief Calculate the rate constant for a set of conditions
    /// @param conditions The current environmental conditions of the chemical system
    /// @return The reaction rate constant
    virtual double Calculate(const Conditions& conditions) const
    {
      return 0.0;
    }

  };

}  // namespace micm