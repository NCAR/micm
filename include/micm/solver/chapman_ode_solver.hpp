/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/solver/solver.hpp>

namespace micm
{

  /**
   * @brief An implementation of the Chapman mechnanism solver
   *
   */
  class ChapmanODESolver : public Solver
  {
   private:
   public:
    /// @brief Default constructor
    ChapmanODESolver() = default;
    /// @brief Move the system to the next state
    /// @param state The collection of species concentrations
    void Solve(double state[]);
  };

}  // namespace micm