/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cstddef>
#include <vector>

namespace micm
{

  /**
   * @brief A base class to represent any time of solver
   *
   */
  class Solver
  {
   protected:
    /// @brief The size of the system
    std::size_t size_;

   public:
    /// @brief Default constructor
    Solver();
    /// @brief Virtual destructor
    virtual ~Solver() = default;

    /// @brief A virtual function to be defined by any solver baseclass
    /// @param time_start Time step to start at
    /// @param time_end Time step to end at
    /// @param number_densities Species concentrations in molecules / cm3
    /// @return A new state representing the species concentrations
    virtual std::vector<double> Solve(double time_start, double time_end, std::vector<double> number_densities) = 0;
  };

  inline Solver::Solver()
      : size_()
  {
  }

}  // namespace micm
