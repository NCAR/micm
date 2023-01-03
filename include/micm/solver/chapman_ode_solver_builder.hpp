/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/solver/solver.hpp>
#include <micm/solver/solver_builder.hpp>

namespace micm
{

  template<typename DataType>
  class ChapmanODESolverBuilder : public SolverBuilder<DataType>
  {
   private:
   public:
    ChapmanODESolverBuilder();

    Solver<DataType> Build();
  };

  template<typename DataType>
  inline ChapmanODESolverBuilder<DataType>::ChapmanODESolverBuilder()
      : SolverBuilder<DataType>()
  {
  }

}  // namespace micm