// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Much of this solver was formulated and implemented from this book:
// Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd
// edition. ed. Springer, Berlin ; New York. The source code for many (all?) of the solvers in that book can be found here:
// http://www.unige.ch/~hairer/software.html
//
// Some extensions to the rosenbrock solver formulated there were formulated in this paper
// Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. Benchmarking stiff ode solvers for
// atmospheric chemistry problems II: Rosenbrock solvers. Atmospheric Environment 31, 3459–3472.
// https://doi.org/10.1016/S1352-2310(97)83212-8
#pragma once

#include <micm/jit/jit_function.hpp>
#include <micm/jit/solver/jit_linear_solver.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/util/random_string.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <chrono>
#include <ctime>
#include <source_location>

namespace micm
{
  struct JitRosenbrockSolverParameters;

  /// @brief A Rosenbrock solver with JIT-compiled optimizations
  template<class RatesPolicy, class LinearSolverPolicy>
  class JitRosenbrockSolver : public AbstractRosenbrockSolver<
                                  RatesPolicy,
                                  LinearSolverPolicy,
                                  JitRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>
  {
    llvm::orc::ResourceTrackerSP function_resource_tracker_;
    using FuncPtr = void (*)(double*, const double);
    FuncPtr alpha_minus_jacobian_ = nullptr;

   public:
    /// @brief Solver parameters typename
    using ParametersType = JitRosenbrockSolverParameters;

    JitRosenbrockSolver(const JitRosenbrockSolver&) = delete;
    JitRosenbrockSolver& operator=(const JitRosenbrockSolver&) = delete;
    JitRosenbrockSolver(JitRosenbrockSolver&& other)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, JitRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(
              std::move(other)),
          function_resource_tracker_(std::move(other.function_resource_tracker_)),
          alpha_minus_jacobian_(std::move(other.alpha_minus_jacobian_))
    {
      other.alpha_minus_jacobian_ = NULL;
    }

    JitRosenbrockSolver& operator=(JitRosenbrockSolver&& other)
    {
      AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, JitRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>::
      operator=(std::move(other));
      function_resource_tracker_ = std::move(other.function_resource_tracker_);
      alpha_minus_jacobian_ = std::move(other.alpha_minus_jacobian_);
      other.alpha_minus_jacobian_ = NULL;
      return *this;
    }

    /// @brief Builds a Rosenbrock solver for the given system and solver parameters
    /// @param parameters Solver parameters
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    /// @param jacobian Jacobian matrix
    JitRosenbrockSolver(
        RosenbrockSolverParameters parameters,
        LinearSolverPolicy linear_solver,
        RatesPolicy rates,
        auto& jacobian)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, JitRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(
              parameters,
              std::move(linear_solver),
              std::move(rates),
              jacobian)
    {
      this->GenerateAlphaMinusJacobian(jacobian);
    }

    ~JitRosenbrockSolver()
    {
      if (function_resource_tracker_ != NULL)
      {
        llvm::ExitOnError exit_on_error;
        exit_on_error(function_resource_tracker_->remove());
      }
    }

    /// @brief compute [alpha * I - dforce_dy]
    /// @param jacobian Jacobian matrix (dforce_dy)
    /// @param alpha
    template<class SparseMatrixPolicy>
    void AlphaMinusJacobian(SparseMatrixPolicy& jacobian, const double& alpha) const
    {
      if (jacobian.GroupVectorSize() != jacobian.NumberOfBlocks())
      {
        std::string msg =
            "JIT functions require the number of grid cells solved together (" +
            std::to_string(jacobian.NumberOfBlocks()) + ") to match the vector dimension template "
            "parameter, currently: " +
            std::to_string(jacobian.GroupVectorSize());
        throw std::system_error(make_error_code(MicmJitErrc::InvalidMatrix), msg);
      }
      double a = alpha;
      if (alpha_minus_jacobian_)
      {
        alpha_minus_jacobian_(jacobian.AsVector().data(), a);
      }
      else
      {
        throw std::system_error(
            make_error_code(MicmJitErrc::MissingJitFunction), std::source_location::current().function_name());
      }
    }

   private:
    void GenerateAlphaMinusJacobian(auto& jacobian)
    {
      auto diagonal_elements = jacobian.DiagonalIndices(0);
      // save sizes needed throughout the function
      std::size_t n_cells = jacobian.GroupVectorSize();
      std::size_t number_of_nonzero_jacobian_elements = jacobian.AsVector().size();

      // Create the JitFunction with the modified name
      JitFunction func = JitFunction::Create()
                             .SetName("alpha_minus_jacobian")
                             .SetArguments({ { "jacobian", JitType::DoublePtr }, { "alpha", JitType::Double } })
                             .SetReturnType(JitType::Void);

      // constants
      llvm::Value* zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
      llvm::Value* negative_one = llvm::ConstantFP::get(*(func.context_), llvm::APFloat(-1.0));
      llvm::Value* jacobian_size =
          llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, number_of_nonzero_jacobian_elements));

      // types
      llvm::Type* double_type = func.GetType(JitType::Double);

      // iterative over the blocks of the jacobian and add the alpha value
      // jacobian_vector[i_elem + i_cell] += alpha;
      for (const auto& i_elem : diagonal_elements)
      {
        llvm::Value* ptr_index[1];

        auto cell_loop = func.StartLoop("add alpha", 0, n_cells);
        llvm::Value* elem_id = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, i_elem));

        ptr_index[0] = func.builder_->CreateNSWAdd(cell_loop.index_, elem_id);

        llvm::Value* indexer =
            func.builder_->CreateGEP(double_type, func.arguments_[0].ptr_, ptr_index, "index jacobian array");
        llvm::Value* jacobian_element = func.builder_->CreateLoad(double_type, indexer, "load jacobian element");

        jacobian_element = func.builder_->CreateFAdd(jacobian_element, func.arguments_[1].ptr_, "add alpha");
        func.builder_->CreateStore(jacobian_element, indexer);

        func.EndLoop(cell_loop);
      }

      func.builder_->CreateRetVoid();

      auto target = func.Generate();

      alpha_minus_jacobian_ = (FuncPtr)(intptr_t)target.second;
      function_resource_tracker_ = target.first;
    }
  };
}  // namespace micm