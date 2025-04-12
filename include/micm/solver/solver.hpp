// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_temporary_variables.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_temporary_variables.hpp>
#include <micm/solver/solver_result.hpp>

#include <type_traits>

namespace micm
{
  template<class SolverPolicy, class StatePolicy>
  class Solver
  {
   private:
    using SolverParametersType = typename SolverPolicy::ParametersType;
    using DenseMatrixType = typename StatePolicy::DenseMatrixPolicyType;

    StateParameters state_parameters_;
    std::vector<micm::Process> processes_;
    System system_;

   public:
    using SolverPolicyType = SolverPolicy;
    using StatePolicyType = StatePolicy;
    SolverPolicy solver_;
    SolverParametersType solver_parameters_;

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        SolverParametersType solver_parameters,
        std::vector<micm::Process> processes,
        System system)
        : solver_(std::move(solver)),
          state_parameters_(state_parameters),
          solver_parameters_(solver_parameters),
          processes_(std::move(processes)),
          system_(std::move(system))
          
    {
    }

    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    Solver(Solver&& other)
        : solver_(std::move(other.solver_)),
          processes_(std::move(other.processes_)),
          state_parameters_(other.state_parameters_),
          solver_parameters_(other.solver_parameters_)
    {
    }
    Solver& operator=(Solver&& other)
    {
      std::swap(this->solver_, other.solver_);
      state_parameters_ = other.state_parameters_;
      solver_parameters_ = other.solver_parameters_;
      std::swap(this->processes_, other.processes_);
      return *this;
    }

    SolverResult Solve(double time_step, StatePolicy& state)
    {
      auto result = solver_.Solve(time_step, state, solver_parameters_);
      state.variables_.Max(0.0);
      return result;
    }

    // Overloaded Solve function to change parameters
    SolverResult Solve(double time_step, StatePolicy& state, const SolverParametersType& params)
    {
      solver_parameters_ = params;
      return solver_.Solve(time_step, state, params);
    }

    /// @brief Returns the number of grid cells
    /// @return
    std::size_t GetNumberOfGridCells() const
    {
      return state_parameters_.number_of_grid_cells_;
    }

    /// @brief Returns the number of species
    /// @return
    std::size_t GetNumberOfSpecies() const
    {
      return state_parameters_.number_of_species_;
    }

    std::size_t GetNumberOfReactions() const
    {
      return state_parameters_.number_of_rate_constants_;
    }

    StatePolicy GetState() const
    {
      auto state = std::move(StatePolicy(state_parameters_));
      if constexpr (std::is_convertible_v<typename SolverPolicy::ParametersType, RosenbrockSolverParameters>)
      {
        state.temporary_variables_ =
            std::make_unique<RosenbrockTemporaryVariables<DenseMatrixType>>(state_parameters_, solver_parameters_);
      }
      else if constexpr (std::is_same_v<typename SolverPolicy::ParametersType, BackwardEulerSolverParameters>)
      {
        state.temporary_variables_ = std::make_unique<BackwardEulerTemporaryVariables<DenseMatrixType>>(state_parameters_);
      }
      else
      {
        throw std::runtime_error(
            "Solver type not supported! Parameter type: " +
            std::string(typeid(typename SolverPolicy::ParametersType).name()));
      }
      return state;
    }

    System GetSystem() const
    {
      return system_;
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
