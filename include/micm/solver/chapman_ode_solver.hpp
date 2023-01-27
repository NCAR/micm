/* Copyright (C) 2023 National Center for Atmospheric Research,
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
    ChapmanODESolver();
    ~ChapmanODESolver();
    /// @brief Move the system to the next state
    /// @param state The collection of species concentrations
    void Solve(double state[]) override;
  };

  inline ChapmanODESolver::ChapmanODESolver()
  {
  }

  inline ChapmanODESolver::~ChapmanODESolver()
  {
  }

  inline void ChapmanODESolver::Solve(double state[])
  {
  }

}  // namespace micm