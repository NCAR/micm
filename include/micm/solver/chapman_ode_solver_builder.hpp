/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/solver.hpp>
#include <micm/solver/solver_builder.hpp>

namespace micm
{

  /**
   * @brief A builder that generates a solver for the Chapman mechanism
   *
   * @tparam T
   */
  template<typename T>
  class ChapmanODESolverBuilder : public SolverBuilder<T>
  {
   private:
   public:
    /// @brief Default constructor
    ChapmanODESolverBuilder();

    /// @brief Adds a micm::Process to this solver
    /// @param process Some process
    /// @return A reference to this solver builder
    SolverBuilder<T>& For(Process process);
    /// @brief Adds zero or more processes to this solver
    /// @param process Some processes a vector of processes
    /// @return A reference to this solver builder
    SolverBuilder<T>& For(std::vector<Process> processes);
    /// @brief Returns the final
    /// @return
    std::unique_ptr<Solver<T>> Build();
  };

  template<typename T>
  inline ChapmanODESolverBuilder<T>::ChapmanODESolverBuilder()
      : SolverBuilder<T>()
  {
  }

  template<typename T>
  inline SolverBuilder<T>& ChapmanODESolverBuilder<T>::For(Process process)
  {
    // TODO: insert return statement here
  }

  template<typename T>
  inline SolverBuilder<T>& ChapmanODESolverBuilder<T>::For(std::vector<Process> processes)
  {
    // TODO: insert return statement here
  }

  template<typename T>
  inline std::unique_ptr<Solver<T>> ChapmanODESolverBuilder<T>::Build()
  {
    return std::make_unique<ChapmanODESolver<T>>();
  }

}  // namespace micm