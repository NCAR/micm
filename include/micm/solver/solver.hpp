// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/process/reaction_rate_store.hpp>
#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_temporary_variables.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_temporary_variables.hpp>
#include <micm/solver/solver_result.hpp>
#include <micm/util/matrix.hpp>

#include <algorithm>
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
    std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixType&)>>
        update_state_parameters_functions_;
    ReactionRateConstantStore store_;
    std::vector<std::function<void(const DenseMatrixType&, DenseMatrixType&)>> initialize_constraint_parameters_functions_;

   public:
    using SolverPolicyType = SolverPolicy;
    using StatePolicyType = StatePolicy;
    SolverPolicy solver_;
    SolverParametersType solver_parameters_;

    Solver(const Solver&) = delete;

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        SolverParametersType solver_parameters,
        std::vector<micm::Process> processes,
        System system)
        : Solver(
              std::move(solver),
              state_parameters,
              solver_parameters,
              std::move(processes),
              std::move(system),
              {})
    {
    }

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        SolverParametersType solver_parameters,
        std::vector<micm::Process> processes,
        System system,
        const std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixType&)>>&
            update_state_parameters_functions)
        : solver_(std::move(solver)),
          state_parameters_(state_parameters),
          solver_parameters_(solver_parameters),
          processes_(std::move(processes)),
          system_(std::move(system)),
          update_state_parameters_functions_(update_state_parameters_functions),
          store_(ReactionRateConstantStore::BuildFrom(processes_))
    {
      if constexpr (requires { solver_.rates_.BuildCudaStore(store_); })
        solver_.rates_.BuildCudaStore(store_);
    }

    Solver(
        SolverPolicy&& solver,
        StateParameters state_parameters,
        SolverParametersType solver_parameters,
        std::vector<micm::Process> processes,
        System system,
        const std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixType&)>>&
            update_state_parameters_functions,
        const std::vector<std::function<void(const DenseMatrixType&, DenseMatrixType&)>>&
            initialize_constraint_parameters_functions)
        : solver_(std::move(solver)),
          state_parameters_(state_parameters),
          solver_parameters_(solver_parameters),
          processes_(std::move(processes)),
          system_(std::move(system)),
          update_state_parameters_functions_(update_state_parameters_functions),
          store_(ReactionRateConstantStore::BuildFrom(processes_)),
          initialize_constraint_parameters_functions_(initialize_constraint_parameters_functions)
    {
      if constexpr (requires { solver_.rates_.BuildCudaStore(store_); })
        solver_.rates_.BuildCudaStore(store_);
    }

    Solver(Solver&& other)
        : solver_(std::move(other.solver_)),
          processes_(std::move(other.processes_)),
          state_parameters_(other.state_parameters_),
          solver_parameters_(other.solver_parameters_),
          system_(std::move(other.system_)),
          update_state_parameters_functions_(std::move(other.update_state_parameters_functions_)),
          store_(std::move(other.store_)),
          initialize_constraint_parameters_functions_(std::move(other.initialize_constraint_parameters_functions_))
    {
    }

    Solver& operator=(const Solver&) = delete;

    Solver& operator=(Solver&& other)
    {
      std::swap(this->solver_, other.solver_);
      state_parameters_ = other.state_parameters_;
      solver_parameters_ = other.solver_parameters_;
      std::swap(this->processes_, other.processes_);
      std::swap(this->system_, other.system_);
      std::swap(this->update_state_parameters_functions_, other.update_state_parameters_functions_);
      std::swap(this->store_, other.store_);
      std::swap(this->initialize_constraint_parameters_functions_, other.initialize_constraint_parameters_functions_);
      return *this;
    }

    SolverResult Solve(double time_step, StatePolicy& state)
    {
      for (const auto& init_func : initialize_constraint_parameters_functions_)
        init_func(state.variables_, state.custom_rate_parameters_);
      auto result = solver_.Solve(time_step, state, solver_parameters_);
      PostSolveClamp(state);
      return result;
    }

    // Overloaded Solve function to change parameters
    SolverResult Solve(double time_step, StatePolicy& state, const SolverParametersType& params)
    {
      solver_parameters_ = params;
      for (const auto& init_func : initialize_constraint_parameters_functions_)
        init_func(state.variables_, state.custom_rate_parameters_);
      auto result = solver_.Solve(time_step, state, params);
      PostSolveClamp(state);
      return result;
    }

    /// @brief Returns the maximum number of grid cells per state
    /// @return Number of grid cells
    /// @details This is the maximum number of grid cells that can fit
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

    StatePolicy GetState(const std::size_t number_of_grid_cells = 1) const
    {
      StatePolicy state(state_parameters_, number_of_grid_cells);

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

    const ReactionRateConstantStore& GetRateConstantStore() const
    {
      return store_;
    }

    /// @brief Update state parameters based on current conditions (temperature, pressure, etc.)
    ///        Invokes all registered parameter update functions for external models and constraints
    ///        to recompute temperature-dependent values (e.g., aerosol rate constants, equilibrium constants)
    ///        Should be called before solving if conditions have changed since the last solve
    /// @param state State object containing conditions and custom_rate_parameters to be updated
    void UpdateStateParameters(StatePolicy& state)
    {
      // External update functions must run first: they populate custom_rate_parameters_
      // which user-defined and surface rate constants read during calculation.
      for (const auto& update_func : update_state_parameters_functions_)
        update_func(state.conditions_, state.custom_rate_parameters_);

      // Dispatch to GPU path if the RatesPolicy (e.g. CudaProcessSet) exposes
      // GpuCalculateRateConstants; otherwise use the CPU path.
      if constexpr (requires { solver_.rates_.GpuCalculateRateConstants(store_, state); })
        solver_.rates_.GpuCalculateRateConstants(store_, state);
      else
      {
        ReactionRateConstantStore::EvaluateCpuRateConstants(store_, state);
        ReactionRateConstantStore::CpuCalculateRateConstants(store_, state);
      }
    }

    LambdaRateConstantParameters& GetLambdaRateConstantByName(const std::string& name)
    {
      return const_cast<LambdaRateConstantParameters&>(std::as_const(*this).GetLambdaRateConstantByName(name));
    }

    const LambdaRateConstantParameters& GetLambdaRateConstantByName(const std::string& name) const
    {
      for (const auto& process : processes_)
      {
        if (const auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
        {
          if (const auto* params = std::get_if<LambdaRateConstantParameters>(&reaction->rate_constant_))
          {
            if (params->label_ == name)
              return *params;
          }
        }
      }
      throw MicmException(
          MicmSeverity::Error,
          MICM_ERROR_CATEGORY_SOLVER,
          MICM_SOLVER_ERROR_CODE_RATE_CONSTANT_NOT_FOUND,
          "Lambda rate constant with name '" + name + "' not found in any process");
    }

   private:
    /// @brief Clamp state variables to non-negative after a solve
    ///        For DAE systems, only ODE variables are clamped; algebraic variables are left unclamped
    void PostSolveClamp(StatePolicy& state)
    {
      if (state.constraint_size_ > 0)
      {
        for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
          for (std::size_t i_var = 0; i_var < state.variables_.NumColumns(); ++i_var)
            if (state.upper_left_identity_diagonal_[i_var] > 0.0)
              state.variables_[i_cell][i_var] = std::max(0.0, state.variables_[i_cell][i_var]);
      }
      else
      {
        state.variables_.Max(0.0);
      }
    }
  };

}  // namespace micm
