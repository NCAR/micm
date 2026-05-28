// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>

namespace micm
{
  /// @brief The final state the solver was in after the Solve function finishes
  enum class SolverState
  {
    /// @brief This is the initial value at the start of the Solve function
    NOT_YET_CALLED,
    /// @brief This is only used for control flow in the Solve function
    RUNNING,
    /// @brief A successful integration will have this value
    CONVERGED,
    /// @brief If the number of steps exceeds the maximum value on the solver parameter, this value will be returned
    CONVERGENCE_EXCEEDED_MAX_STEPS,
    /// @brief Very stiff systems will likely result in a step size that is not useable for the solver
    STEP_SIZE_TOO_SMALL,
    /// @brief Matrices that are singular more than once will set this value. At present, this should never be returned
    REPEATEDLY_SINGULAR_MATRIX,
    /// @brief Mostly this value is returned by systems that tend toward chemical explosions
    NA_N_DETECTED,
    /// @brief Can happen when unititialized memory is used in the solver
    INF_DETECTED,
    /// @brief Used for backward euler. This allows us to "succeed" in the same way that cam-chem does
    ACCEPTING_UNCONVERGED_INTEGRATION,
    /// @brief Newton iteration to initialize algebraic constraint variables failed to converge
    CONSTRAINT_INITIALIZATION_FAILED
  };

  struct SolverStats
  {
    /// @brief The number of forcing function calls
    uint64_t function_calls_{};
    /// @brief The number of jacobian function calls
    uint64_t jacobian_updates_{};
    /// @brief The total number of internal time steps taken
    uint64_t number_of_steps_{};
    /// @brief The number of accepted integrations
    uint64_t accepted_{};
    /// @brief The number of rejected integrations
    uint64_t rejected_{};
    /// @brief The number of LU decompositions
    uint64_t decompositions_{};
    /// @brief The number of linear solves
    uint64_t solves_{};
    /// @brief The number of constraint initialization iterations performed
    uint64_t constraint_init_iterations_{};
    /// @brief The final time the solver iterated to
    double final_time_{};
  };

  inline std::string SolverStateToString(const SolverState& state)
  {
    switch (state)
    {
      case SolverState::NOT_YET_CALLED: return "Not Yet Called";
      case SolverState::RUNNING: return "Running";
      case SolverState::CONVERGED: return "Converged";
      case SolverState::CONVERGENCE_EXCEEDED_MAX_STEPS: return "Convergence Exceeded Max Steps";
      case SolverState::STEP_SIZE_TOO_SMALL: return "Step Size Too Small";
      case SolverState::REPEATEDLY_SINGULAR_MATRIX: return "Repeatedly Singular Matrix";
      case SolverState::NA_N_DETECTED: return "NaNDetected";
      case SolverState::INF_DETECTED: return "InfDetected";
      case SolverState::ACCEPTING_UNCONVERGED_INTEGRATION: return "AcceptingUnconvergedIntegration";
      case SolverState::CONSTRAINT_INITIALIZATION_FAILED: return "Constraint Initialization Failed";
      default: return "Unknown";
    }
  }

  struct [[nodiscard]] SolverResult
  {
    /// @brief The final state the solver was in
    SolverState state_ = SolverState::NOT_YET_CALLED;
    /// @brief A collection of runtime state for this call of the solver
    SolverStats stats_{};
  };
}  // namespace micm
