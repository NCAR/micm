/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/solver/backward_euler.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

#include <memory>
#include <variant>

namespace micm
{
  using SolverParameters = std::variant<std::monostate, micm::RosenbrockSolverParameters, micm::BackwardEulerSolverParameters>;

  class SolverImplBase
  {
   public:
    virtual ~SolverImplBase() = default;
  };

  template<class LinearSolverPolicy, class ProcessSetPolicy>
  struct SolverImpl : public SolverImplBase
  {
    ProcessSetPolicy process_set_;
    LinearSolverPolicy linear_solver_;
  };

  class Solver
  {
   private:
    SolverParameters parameters_;
    std::size_t number_of_grid_cells_;
    std::size_t number_of_species_;
    std::size_t number_of_reactions_;
    std::unique_ptr<SolverImplBase> solver_impl_;

   public:
    Solver()
        : parameters_(std::monostate{}),
          number_of_grid_cells_(0),
          number_of_species_(0),
          number_of_reactions_(0)
    {
    }

    Solver(
        SolverImplBase* solver_impl,
        SolverParameters parameters,
        std::size_t number_of_grid_cells,
        std::size_t number_of_species,
        std::size_t number_of_reactions)
        : solver_impl_(solver_impl),
          parameters_(parameters),
          number_of_grid_cells_(number_of_grid_cells),
          number_of_species_(number_of_species),
          number_of_reactions_(number_of_reactions)
    {
    }

    void Solve()
    {
      if (std::holds_alternative<micm::RosenbrockSolverParameters>(parameters_))
      {
        // call Rosenbrock solver
      }
      else if (std::holds_alternative<micm::BackwardEulerSolverParameters>(parameters_))
      {
        micm::BackwardEuler be;
        // be.solve()
      }
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
  };
}  // namespace micm
