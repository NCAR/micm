/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/solver.hpp>
#include <micm/system/system.hpp>
#include <vector>

namespace micm
{

  template<typename DataType>
  class SolverBuilder
  {
   protected:
    const System<DataType> system_;
    const std::vector<Process> processes_;

   public:
    SolverBuilder();
    SolverBuilder(System<DataType> system);

    SolverBuilder& For(Process process);
    SolverBuilder& For(std::vector<Process> processes);

    virtual Solver<DataType> Build() = 0;
  };

  template<typename DataType>
  inline SolverBuilder<DataType>::SolverBuilder()
    : system_()
  {
  }
}  // namespace micm