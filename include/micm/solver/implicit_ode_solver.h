/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_IMPLICIT_ODE_SOLVER_H
#define MICM_IMPLICIT_ODE_SOLVER_H

#include <micm/solver/implicit_ode_solver_builder.h>
#include <micm/solver/solver.h>

template <typename T> class ImplicitODESolver : public Solver<T> {
  private:

  public:
    ImplicitODESolver(ImplicitODESolverBuilder<T>);

    void Solve(T[]);
};

#endif