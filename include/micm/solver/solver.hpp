// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/solver_result.hpp>

namespace micm
{
  template<class SolverPolicy, class StatePolicy>
  class Solver
  {
   private:
    std::size_t number_of_grid_cells_;
    std::size_t number_of_species_;
    std::size_t number_of_reactions_;
    StateParameters state_parameters_;
    std::vector<micm::Process> processes_;

   public:
    SolverPolicy solver_;

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        std::size_t number_of_grid_cells,
        std::size_t number_of_species,
        std::size_t number_of_reactions,
        std::vector<micm::Process> processes)
        : solver_(std::move(solver)),
          number_of_grid_cells_(number_of_grid_cells),
          number_of_species_(number_of_species),
          number_of_reactions_(number_of_reactions),
          state_parameters_(state_parameters),
          processes_(std::move(processes))
    {
    }

    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    Solver(Solver&& other)
        : solver_(std::move(other.solver_)),
          processes_(std::move(other.processes_)),
          number_of_grid_cells_(other.number_of_grid_cells_),
          number_of_species_(other.number_of_species_),
          number_of_reactions_(other.number_of_reactions_),
          state_parameters_(other.state_parameters_)
    {
    }
    Solver& operator=(Solver&& other)
    {
      std::swap(this->solver_, other.solver_);
      number_of_grid_cells_ = other.number_of_grid_cells_;
      number_of_species_ = other.number_of_species_;
      number_of_reactions_ = other.number_of_reactions_;
      state_parameters_ = other.state_parameters_;
      std::swap(this->processes_, other.processes_);
      return *this;
    }

    SolverResult Solve(double time_step, StatePolicy& state)
    {
      return solver_.Solve(time_step, state);
    }

    /// @brief Returns the number of grid cells
    /// @return
    std::size_t GetNumberOfGridCells() const
    {
      return number_of_grid_cells_;
    }

    /// @brief Returns the number of species
    /// @return
    std::size_t GetNumberOfSpecies() const
    {
      return number_of_species_;
    }

    /// @brief Returns the number of reactions
    /// @return
    std::size_t GetNumberOfReactions() const
    {
      return number_of_reactions_;
    }

    StatePolicy GetState() const
    {
      return StatePolicy(state_parameters_);
    }

    std::vector<micm::Process> GetProcesses() const
    {
      return processes_;
    }

    void CalculateRateConstants(StatePolicy& state)
    {
      Process::CalculateRateConstants(processes_, state);
    }
  };

}  // namespace micm
