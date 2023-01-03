/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/solver/solver.hpp>

namespace micm
{

  template<typename DataType>
  class ChapmanODESolver : public Solver<DataType>
  {
   private:
   public:
    void Solve(DataType);
  };

}  // namespace micm