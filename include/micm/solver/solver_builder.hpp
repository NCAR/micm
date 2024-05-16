/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/backward_euler_solver_parameters.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/solver.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <system_error>

namespace micm
{
  class SolverBuilder
  {
   protected:
    micm::System system_;
    std::size_t number_of_grid_cells_;
    std::vector<micm::Process> reactions_;
    SolverParameters options_;
    bool ignore_unused_species_ = true;
    bool reorder_state_ = true;

   public:
    /// @brief Set the chemical system
    /// @param system
    /// @return
    SolverBuilder& SetSystem(micm::System system);

    /// @brief Set the reactions
    /// @param reactions
    /// @return
    SolverBuilder& SetReactions(std::vector<micm::Process> reactions);

    /// @brief Set the number of grid cells
    /// @param number_of_grid_cells
    /// @return
    SolverBuilder& SetNumberOfGridCells(int number_of_grid_cells);

    /// @brief Choose a rosenbrock solver
    /// @param options
    /// @return
    SolverBuilder& SetSolverParameters(const RosenbrockSolverParameters& options);

    /// @brief Choose a backward euler solver
    /// @param options
    /// @return
    SolverBuilder& SetSolverParameters(const BackwardEulerSolverParameters& options);

    /// @brief
    /// @return
    template<class MatrixPolicy, class SparseMatrixPolicy>
    Solver<State<MatrixPolicy, SparseMatrixPolicy>> Build();

   protected:
    /// @brief
    /// @return
    virtual Solver BuildBackwardEulerSolver() = 0;
    virtual ~SolverBuilder() = default;

    /// @brief
    /// @tparam ProcessSetPolicy
    /// @return
    template<class ProcessSetPolicy>
    void UnusedSpeciesCheck();

    /// @brief Get a species map properly ordered
    /// @return
    template<class MatrixPolicy, class ProcessSetPolicy>
    std::map<std::string, std::size_t> GetSpeciesMap() const;

    /// @brief Set the absolute tolerances per species
    /// @param parameters
    /// @param species_map
    /// @return
    void SetAbsoluteTolerances(std::vector<double>& tolerances, const std::map<std::string, std::size_t>& species_map) const;

    /// @brief Return the labels of the custom parameters
    /// @return
    std::vector<std::string> GetCustomParameterLabels() const;

    std::vector<std::size_t> GetJacobianDiagonalElements(auto jacobian) const;
  };

  template<class MatrixPolicy, class SparseMatrixPolicy>
  class CpuSolverBuilder : public SolverBuilder
  {
   public:
    Solver BuildBackwardEulerSolver() override;
  };

  template<class MatrixPolicy, class SparseMatrixPolicy>
  class GpuSolverBuilder : public SolverBuilder
  {
   public:
    Solver BuildBackwardEulerSolver() override;
  };

}  // namespace micm

#include "solver_builder.inl"