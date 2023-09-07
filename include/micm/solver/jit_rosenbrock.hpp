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
    using FuncPtr = void (*)(const double*, const double);
    FuncPtr alpha_minus_jacobian_ = nullptr;

   public:
    /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
    /// @param system The chemical system to create the solver for
    /// @param processes The collection of chemical processes that will be applied during solving
    JitRosenbrockSolver(
        std::shared_ptr<JitCompiler> compiler,
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters)
        : RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>(system, processes, parameters),
          compiler_(compiler)
    {
      this->GenerateAlphaMinusJacobian();
    }

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
    {
      if (alpha_minus_jacobian_)
      {
        for(auto& elem: jacobian.AsVector()) {
          std::cout << elem << ", ";
        }
        std::cout << std::endl;
        alpha_minus_jacobian_(jacobian.AsVector().data(), alpha);
        for(auto& elem: jacobian.AsVector()) {
          std::cout << elem << ", ";
        }
        std::cout << std::endl;
      }
    }

   private:
    /*
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
    */
    void GenerateAlphaMinusJacobian()
    {
      // save sizes needed throughout the function
      std::size_t L = this->jacobian_.GroupVectorSize();
      std::size_t number_of_nonzero_jacobian_elements = this->jacobian_.AsVector().size();

      JitFunction func = JitFunction::create(compiler_)
                             .name("alpha_minus_jacobian")
                             .arguments({ { "jacobian", JitType::DoublePtr }, { "alpha", JitType::Double } })
                             .return_type(JitType::Void);

      // constants
      llvm::Value* zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
      llvm::Value* negative_one = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(-1.0));
      llvm::Value* jacobian_size = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, number_of_nonzero_jacobian_elements));

      // types
      llvm::Type* double_type = func.GetType(JitType::Double);

      // multiply each jacobian element by negative
      {
        llvm::Value* ptr_index[1];

        llvm::Value* index_list[2];
        index_list[0] = zero;
        auto loop = func.StartLoop("negate jacobian", 0, number_of_nonzero_jacobian_elements);
        index_list[1] = loop.index_;
        ptr_index[0] = loop.index_;
        llvm::Value* indexer = func.builder_->CreateInBoundsGEP(double_type, func.arguments_[0].ptr_, ptr_index, "create an indexer into the jacobian array");
        llvm::Value* jacobian_element = func.builder_->CreateLoad(double_type, indexer, "load jacobian element");
        jacobian_element = func.builder_->CreateFMul(jacobian_element, negative_one, "negate");
        func.builder_->CreateStore(jacobian_element, indexer);

        func.EndLoop(loop);
      }

      func.builder_->CreateRetVoid();

      func.module_->print(llvm::outs(), nullptr);
      auto target = func.Generate();

      alpha_minus_jacobian_ = (FuncPtr)(intptr_t)target.second;
      function_resource_tracker_ = target.first;
    }
  };
}  // namespace micm