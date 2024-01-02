/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cstddef>
#include <iterator>
#include <memory>
#include <micm/system/conditions.hpp>
#include <vector>

namespace micm
{
  /// @brief A base class for any type of rate constant
  class RateConstant
  {
   public:
    /// @brief Virtual destructor
    virtual ~RateConstant(){};

    /// @brief Deep copy
    virtual std::unique_ptr<RateConstant> clone() const = 0;

    /// @brief Returns a set of labels for user-defined rate constant parameters
    /// @return Vector of custom parameter labels
    virtual std::vector<std::string> CustomParameters() const
    {
      return std::vector<std::string>{};
    }

    /// @brief Returns the number of custom parameters
    /// @return Number of custom parameters
    virtual std::size_t SizeCustomParameters() const
    {
      return 0;
    }

    /// @brief Calculate the rate constant for a set of conditions
    /// @param conditions The current environmental conditions of the chemical system
    /// @return The reaction rate constant
    virtual double calculate(const Conditions& conditions) const
    {
      return 0;
    }

    /// @brief Calculate the rate constant for a set of conditions
    /// @param conditions The current environmental conditions of the chemical system
    /// @param custom_parameters User defined rate constant parameters
    /// @return The reaction rate constant
    virtual double calculate(const Conditions& conditions, std::vector<double>::const_iterator custom_parameters) const
    {
      return 0;
    }
  };

}  // namespace micm