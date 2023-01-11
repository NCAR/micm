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
   */
  class ChapmanODESolverBuilder : public SolverBuilder
  {
   private:
   public:
    /// @brief Default constructor
    ChapmanODESolverBuilder();

    /// @brief Adds a micm::Process to this solver
    /// @param process Some process
    /// @return A reference to this solver builder
    SolverBuilder& For(Process process);
    /// @brief Adds zero or more processes to this solver
    /// @param process Some processes a vector of processes
    /// @return A reference to this solver builder
    SolverBuilder& For(std::vector<Process> processes);
    /// @brief Returns the final
    /// @return
    std::unique_ptr<Solver> Build();
  };

  inline ChapmanODESolverBuilder::ChapmanODESolverBuilder()
      : SolverBuilder()
  {
  }

  inline SolverBuilder& ChapmanODESolverBuilder::For(Process process)
  {
    // TODO: insert return statement here
  }

  inline SolverBuilder& ChapmanODESolverBuilder::For(std::vector<Process> processes)
  {
    // TODO: insert return statement here
  }

  inline std::unique_ptr<Solver> ChapmanODESolverBuilder::Build()
  {
    return std::make_unique<ChapmanODESolver>();
  }

}  // namespace micm