/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <cstddef>
#include <vector>

namespace micm
{

  /**
   * @brief A base class to represent any time of solver
   *
   */
  class Solver
  {
   protected:
    /// @brief The size of the system
    std::size_t size_;

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
      uint64_t forcing_function_calls{};  // Nfun
      uint64_t jacobian_updates{};        // Njac
      uint64_t number_of_steps{};         // Nstp
      uint64_t accepted{};                // Nacc
      uint64_t rejected{};                // Nrej
      uint64_t decompositions{};          // Ndec
      uint64_t solves{};                  // Nsol
      uint64_t singular{};                // Nsng
      uint64_t total_steps{};             // Ntotstp

      void reset()
      {
        forcing_function_calls = 0;
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

    struct [[nodiscard]] SolverResult
    {
      std::vector<double> result_{};
      SolverState state_ = SolverState::NotYetCalled;
      Rosenbrock_stats stats_{};
      double T{};
    };

   public:
    /// @brief Default constructor
    Solver();
    /// @brief Virtual destructor
    virtual ~Solver() = default;

    /// @brief A virtual function to be defined by any solver baseclass
    /// @param time_start Time step to start at
    /// @param time_end Time step to end at
    /// @param number_densities Species concentrations in molecules / cm3
    /// @param number_density_air The number density of air in molecules / cm3
    /// @return A struct containing results and a status code
    virtual SolverResult Solve(
        const double& time_start,
        const double& time_end,
        const std::vector<double>& number_densities,
        const double& number_density_air) = 0;
  };

  inline Solver::Solver()
      : size_()
  {
  }

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
