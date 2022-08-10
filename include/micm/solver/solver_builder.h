/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SOLVER_BUILDER_H
#define MICM_SOLVER_BUILDER_H

#include <micm/solver/solver.h>
#include <micm/system/system.h>
#include <micm/process/process.h>

#include <vector>

template <typename T> class SolverBuilder {
  private:
    const System system_;
    const std::vector<Process> processes_;
    const std::vector<Solver<T>> embedded_solvers_;

  public:
    SolverBuilder(System system);

    SolverBuilder For(Process process);
    SolverBuilder For(std::vector<Prcoess> processes);
    SolverBuilder Embed(Solver<T> solver);
    virtual Solver<T> Build() = 0;
};

#endif