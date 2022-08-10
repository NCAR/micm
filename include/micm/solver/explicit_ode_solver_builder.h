/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_EXPLICIT_ODE_SOLVER_BUILDER_H
#define MICM_EXPLICIT_ODE_SOLVER_BUILDER_H

#include <micm/solver/solver.h>
#include <micm/solver/solver_builder.h>

template <typename T> class ExplicitODESolverBuilder : public SolverBuilder<T> {
  private:

  public:
    Solver<T> Build();
};

#endif