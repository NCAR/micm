// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_temporary_variables.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_temporary_variables.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/util/matrix.hpp>

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
          solver_parameters_(other.solver_parameters_),
          system_(std::move(other.system_))
    {
    }
    Solver& operator=(Solver&& other)
    {
      std::swap(this->solver_, other.solver_);
      state_parameters_ = other.state_parameters_;
      solver_parameters_ = other.solver_parameters_;
      std::swap(this->processes_, other.processes_);
      std::swap(this->system_, other.system_);
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

    /// @brief Returns the maximum number of grid cells per state
    /// @return Number of grid cells
    /// @details This is the maximum number of grid cells that can that fit
    ///          within one group for vectorized solvers. For non-vectorized solvers,
    ///          there is no limit other than the maximum size of a std::size_t.
    std::size_t MaximumNumberOfGridCells() const
    {
      if constexpr (VectorizableDense<DenseMatrixType>)
      {
        return DenseMatrixType::GroupVectorSize();
      }
      else
      {
        return std::numeric_limits<std::size_t>::max();
      }
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

    StatePolicy GetState(const std::size_t number_of_grid_cells = 1) const
    {
      auto state = std::move(StatePolicy(state_parameters_, number_of_grid_cells));
      if constexpr (std::is_convertible_v<typename SolverPolicy::ParametersType, RosenbrockSolverParameters>)
      {
        state.temporary_variables_ = std::make_unique<RosenbrockTemporaryVariables<DenseMatrixType>>(
            state_parameters_, solver_parameters_, number_of_grid_cells);
      }
      else if constexpr (std::is_same_v<typename SolverPolicy::ParametersType, BackwardEulerSolverParameters>)
      {
        state.temporary_variables_ =
            std::make_unique<BackwardEulerTemporaryVariables<DenseMatrixType>>(state_parameters_, number_of_grid_cells);
      }
      else
      {
        throw std::runtime_error(
            "Solver type not supported! Parameter type: " +
            std::string(typeid(typename SolverPolicy::ParametersType).name()));
      }
      return state;
    }

    const System& GetSystem()
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
