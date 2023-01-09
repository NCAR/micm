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
   * @tparam T The underlying data format of the system
   */
  template<typename T>
  class Solver
  {
   protected:
    /// @brief The size of the system
    std::size_t size_;

   public:
    /// @brief Default constructor
    Solver();
    /// @brief A virtual function to be defined by any solver baseclass
    /// @param T The current species concentrations of the system
    virtual void Solve(T) = 0;
  };

  template<typename T>
  inline Solver<T>::Solver()
    : size_()
  {
  }

}  // namespace micm
