/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_EXPLICIT_ODE_SOLVER_H
#define MICM_EXPLICIT_ODE_SOLVER_H

#include <micm/solver/explicit_ode_solver_builder.h>
#include <micm/solver/solver.h>

template <typename T> class ExplicitODESolver : public Solver<T> {
  private:

  public:
    ExplicitODESolver(ExplicitODESolverBuilder<T>);

    void Solve(T[]);
};

#endif