// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/solver/backward_euler_temporary_variables.hpp>
#include <micm/solver/rosenbrock_temporary_variables.hpp>
#include <variant>

namespace micm
{
  template<class SolverPolicy, class StatePolicy>
  class Solver
  {
   private:
    using SolverParametersType = typename SolverPolicy::ParametersType;
    using DenseMatrixType = typename StatePolicy::DenseMatrixPolicyType;

    std::size_t number_of_grid_cells_;
    std::size_t number_of_species_;
    std::size_t number_of_reactions_;
    StateParameters state_parameters_;
    SolverParametersType solver_parameters_;
    std::vector<micm::Process> processes_;

   public:
    SolverPolicy solver_;

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        SolverParametersType solver_parameters,
        std::size_t number_of_grid_cells,
        std::size_t number_of_species,
        std::size_t number_of_reactions,
        std::vector<micm::Process> processes)
        : solver_(std::move(solver)),
          number_of_grid_cells_(number_of_grid_cells),
          number_of_species_(number_of_species),
          number_of_reactions_(number_of_reactions),
          state_parameters_(state_parameters),
          solver_parameters_(solver_parameters),
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
          state_parameters_(other.state_parameters_),
          solver_parameters_(other.solver_parameters_)
    {
    }
    Solver& operator=(Solver&& other)
    {
      std::swap(this->solver_, other.solver_);
      number_of_grid_cells_ = other.number_of_grid_cells_;
      number_of_species_ = other.number_of_species_;
      number_of_reactions_ = other.number_of_reactions_;
      state_parameters_ = other.state_parameters_;
      solver_parameters_ = other.solver_parameters_;
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
      auto state = StatePolicy(state_parameters_);
      if (std::is_same<typename SolverPolicy::ParametersType, RosenbrockSolverParameters>::value)
      {
        state.temporary_variables_ = std::make_unique<RosenbrockTemporaryVariables<DenseMatrixType>>(state_parameters_, solver_parameters_);
      }
      else if (std::is_same<typename SolverPolicy::ParametersType, BackwardEulerSolverParameters>::value)
      {
        state.temporary_variables_ = std::make_unique<BackwardEulerTemporaryVariables<DenseMatrixType>>(state_parameters_);
      }
      else
      {
        throw std::runtime_error("Solver type not supported!");
      }
      return state;
    }

    std::vector<micm::Process> GetProcesses() const
    {
      return processes_;
    }

    void CalculateRateConstants(StatePolicy& state)
    {
      Process::CalculateRateConstants<DenseMatrixType>(processes_, state);
    }
  };

}  // namespace micm
