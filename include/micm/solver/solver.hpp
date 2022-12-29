/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <stddef.h>

namespace micm
{

  template<typename T>
  class Solver
  {
   private:
    const size_t size_;

   public:
    virtual void Solve(T) = 0;
  };

}  // namespace micm
