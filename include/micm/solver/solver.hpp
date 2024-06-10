// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

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
    SolverPolicy solver_;

   public:
    Solver(
        SolverPolicy solver,
        StateParameters state_parameters,
        std::size_t number_of_grid_cells,
        std::size_t number_of_species,
        std::size_t number_of_reactions)
        : solver_(solver),
          number_of_grid_cells_(number_of_grid_cells),
          number_of_species_(number_of_species),
          number_of_reactions_(number_of_reactions),
          state_parameters_(state_parameters)
    {
    }

    void Solve(double time_step, StatePolicy& state)
    {
      solver_.Solve(time_step, state);
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
  };
}  // namespace micm
