/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cstddef>

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
    virtual void Solve(double state[]) = 0;
  };

  inline Solver::Solver()
      : size_()
  {
  }

}  // namespace micm
