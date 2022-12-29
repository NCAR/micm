/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

namespace micm
{

#include <micm/solver/solver.hpp>
#include <micm/system/system.hpp>
#include <vector>

  template<typename T>
  class SolverBuilder
  {
   private:
    const System system_;
    // const std::vector<Process> processes_;
    const std::vector<Solver<T>> embedded_solvers_;

   public:
    SolverBuilder(System system);

    SolverBuilder& For(Process process);
    SolverBuilder& For(std::vector<Process> processes);
    SolverBuilder& Embed(Solver<T> solver);
    virtual Solver<T> Build() = 0;
  };

  template<typename T>
  inline SolverBuilder<T>::SolverBuilder(System system)
  {
  }

  template<typename T>
  inline SolverBuilder& SolverBuilder<T>::Embed(Solver<T> solver)
  {
    return *this;
  }
}  // namespace micm