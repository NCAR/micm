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
    /// @param state The current species concentrations of the system
    virtual std::vector<double> Solve(std::vector<double> LU, std::vector<double> b) = 0;
  };

  inline Solver::Solver()
      : size_()
  {
  }

}  // namespace micm
