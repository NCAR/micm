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
   * @tparam T The underlying data type of the system
   */
  template<typename T>
  class ChapmanODESolver : public Solver<T>
  {
   private:
   public:
    /// @brief Default constructor
    ChapmanODESolver();
    /// @brief Move the system to the next state
    /// @param The collection of species concentrations
    void Solve(T);
  };

  template<typename T>
  inline ChapmanODESolver<T>::ChapmanODESolver()
    : Solver<T>()
  {
  }

  template<typename T>
  inline void ChapmanODESolver<T>::Solve(T)
  {
  }

}  // namespace micm