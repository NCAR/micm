/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Much of this solver was formulated and implemented from this book:
 * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
 * edition. ed. Springer, Berlin ; New York. The source code for many (all?) of the solvers in that book can be found here:
 * http://www.unige.ch/~hairer/software.html
 *
 * Some extensions to the rosenbrock solver formulated there were formulated in this paper
 * Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. Benchmarking stiff ode solvers for
 * atmospheric chemistry problems II: Rosenbrock solvers. Atmospheric Environment 31, 3459–3472.
 * https://doi.org/10.1016/S1352-2310(97)83212-8
 *
 */
#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/rosenbrock.hpp>

namespace micm
{

  template<template<class> class MatrixPolicy = Matrix, template<class> class SparseMatrixPolicy = SparseMatrix>
  class JitRosenbrockSolver : public RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP function_resource_tracker_;
    void (*alpha_minus_jacobian)(const double *, const double *, double *);

   public:
    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    JitRosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters) 
          : RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>(system, processes, parameters)
        {

        }

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const 
    {
    }

   private:
    void GenerateAlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const 
    {
      const std::size_t n_cells = jacobian.GroupVectorSize();
      for (auto& elem : jacobian.AsVector())
        elem = -elem;
      for (std::size_t i_group = 0; i_group < jacobian.NumberOfGroups(jacobian.size()); ++i_group)
      {
        auto jacobian_vector =
            std::next(jacobian.AsVector().begin(), i_group * jacobian.GroupSize(jacobian.FlatBlockSize()));
        for (const auto& i_elem : this->jacobian_diagonal_elements_)
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            jacobian_vector[i_elem + i_cell] += alpha;
      }
    }
  };
}  // namespace micm