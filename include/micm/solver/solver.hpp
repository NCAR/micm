/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief A base class to represent any time of solver
   *
   */
  class Solver
  {
   public:
    enum class SolverState
    {
      NotYetCalled,
      Converged,
      ConvergenceExceededMaxSteps,
      StepSizeTooSmall,
      RepeatedlySingularMatrix
    };

    struct Rosenbrock_stats
    {
      uint64_t function_calls{};    // Nfun
      uint64_t jacobian_updates{};  // Njac
      uint64_t number_of_steps{};   // Nstp
      uint64_t accepted{};          // Nacc
      uint64_t rejected{};          // Nrej
      uint64_t decompositions{};    // Ndec
      uint64_t solves{};            // Nsol
      uint64_t singular{};          // Nsng
      uint64_t total_steps{};       // Ntotstp

      void reset()
      {
        function_calls = 0;
        jacobian_updates = 0;
        number_of_steps = 0;
        accepted = 0;
        rejected = 0;
        decompositions = 0;
        solves = 0;
        singular = 0;
        total_steps = 0;
      }
    };

    template<typename R>
    struct [[nodiscard]] SolverResult
    {
      /// @brief The new state computed by the solver
      R result_{};
      /// @brief The finals state the solver was in
      SolverState state_ = SolverState::NotYetCalled;
      /// @brief A collection of runtime state for this call of the solver
      Rosenbrock_stats stats_{};
      /// @brief The final time the solver iterated to
      double T{};
    };

   public:
  };

  std::string state_to_string(const Solver::SolverState& state)
  {
    switch (state)
    {
      case micm::Solver::SolverState::NotYetCalled: return "Not Yet Called";
      case micm::Solver::SolverState::Converged: return "Converged";
      case micm::Solver::SolverState::ConvergenceExceededMaxSteps: return "Convergence Exceeded Max Steps";
      case micm::Solver::SolverState::StepSizeTooSmall: return "Step Size Too Small";
      case micm::Solver::SolverState::RepeatedlySingularMatrix: return "Repeatedly Singular Matrix";
      default: return "Unknown";
    }
    return "";
  }

}  // namespace micm
