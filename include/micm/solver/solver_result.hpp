/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

namespace micm
{
  /// @brief The final state the solver was in after the Solve function finishes
  enum class SolverState
  {
    /// @brief This is the initial value at the start of the Solve function
    NotYetCalled,
    /// @brief This is only used for control flow in the Solve function
    Running,
    /// @brief A successful integration will have this value
    Converged,
    /// @brief If the number of steps exceeds the maximum value on the solver parameter, this value will be returned
    ConvergenceExceededMaxSteps,
    /// @brief Very stiff systems will likely result in a step size that is not useable for the solver
    StepSizeTooSmall,
    /// @brief Matrices that are singular more than once will set this value. At present, this should never be returned
    RepeatedlySingularMatrix,
    /// @brief Mostly this value is returned by systems that tend toward chemical explosions
    NaNDetected,
    /// @brief Can happen when unititialized memory is used in the solver
    InfDetected,
    /// @brief Used for backward euler. This allows us to "succeed" in the same way that cam-chem does
    AcceptingUnconvergedIntegration
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

    double solve_timing_{};
    double addforcing_timing_{};
    double jacobian_timing_{};
    double lu_decomp_timing_{};
    double lu_solver_timing_{};
    double convergence_check_timing_{};
    
    /// @brief Set all member variables to zero
    void Reset();
  };

  inline void SolverStats::Reset()
  {
    function_calls_ = 0;
    jacobian_updates_ = 0;
    number_of_steps_ = 0;
    accepted_ = 0;
    rejected_ = 0;
    decompositions_ = 0;
    solves_ = 0;
  }

  inline std::string SolverStateToString(const SolverState& state)
  {
    switch (state)
    {
      case SolverState::NotYetCalled: return "Not Yet Called";
      case SolverState::Running: return "Running";
      case SolverState::Converged: return "Converged";
      case SolverState::ConvergenceExceededMaxSteps: return "Convergence Exceeded Max Steps";
      case SolverState::StepSizeTooSmall: return "Step Size Too Small";
      case SolverState::RepeatedlySingularMatrix: return "Repeatedly Singular Matrix";
      case SolverState::NaNDetected: return "NaNDetected";
      case SolverState::InfDetected: return "InfDetected";
      case SolverState::AcceptingUnconvergedIntegration: return "AcceptingUnconvergedIntegration";
      default: return "Unknown";
    }
  }

  struct [[nodiscard]] SolverResult
  {
    /// @brief The final state the solver was in
    SolverState state_ = SolverState::NotYetCalled;
    /// @brief A collection of runtime state for this call of the solver
    SolverStats stats_{};
    /// @brief The final time the solver iterated to
    double final_time_{};
  };
}  // namespace micm
