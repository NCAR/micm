/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/solver.hpp>
#include <micm/system/system.hpp>
#include <vector>
#include <memory>

namespace micm
{

  /**
   * @brief A base class which implemenets the builder patterns to make different types of solvers
   *
   * @tparam T The underlying datatype of the system
   */
  template<typename T>
  class SolverBuilder
  {
   protected:
    /// @brief The system that the builder will use to generate the solver
    const System<T> system_;
    /// @brief //TODO
    const std::vector<Process> processes_;

   public:
    /// @brief Default constuctor
    SolverBuilder();

    /// @brief A virtual function that adds a micm::Process to a solver
    /// @param process Some process
    /// @return A reference to this solver builder
    virtual SolverBuilder& For(Process process) = 0;
    /// @brief A virtual function that adds zero or more processes to a solver
    /// @param process Some processes a vector of processes
    /// @return A reference to this solver builder
    virtual SolverBuilder& For(std::vector<Process> processes) = 0;

    /// @brief A virtual function that returns the final solver based off of the processes added to this builder
    /// @return A concrete implementation of some solver
    virtual std::unique_ptr<Solver<T>> Build() = 0;
  };

  template<typename T>
  inline SolverBuilder<T>::SolverBuilder()
      : system_()
  {
  }
}  // namespace micm