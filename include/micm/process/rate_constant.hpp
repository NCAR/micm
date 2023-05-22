/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once
#include <cstddef>
#include <iterator>
#include <memory>
#include <micm/solver/state.hpp>
#include <vector>

namespace micm
{

  class System;

  /**
   * @brief A base class for any type of rate constant
   *
   */
  class RateConstant
  {
   public:
    /// @brief Virtual destructor
    virtual ~RateConstant(){};
    /// @brief Deep copy
    virtual std::unique_ptr<RateConstant> clone() const = 0;
    /// @brief Returns the number of doubles needed to hold user-defined rate constant parameters
    /// @return Number of user-defined rate constant parameters
    virtual std::size_t SizeCustomParameters() const;
    /// @brief Calculate the rate constant for a set of conditions
    /// @param state The current state of the chemical system
    /// @param custom_parameters User defined rate constant parameters
    /// @return The reaction rate constant
    virtual double calculate(const State& state, std::vector<double>::const_iterator custom_parameters) const
    {
      return 0;
    };

   private:
  };

  inline std::size_t RateConstant::SizeCustomParameters() const
  {
    return 0;
  }

}  // namespace micm