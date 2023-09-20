// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief A general-use block-diagonal sparse-matrix linear solver with JIT-compiled optimizations
  ///
  /// See LinearSolver class description for algorithm details
  /// The template parameter is the number of blocks (i.e. grid cells) in the block-diagonal matrix
  template<std::size_t L, template<class> class SparseMatrixPolicy>
  class JitLinearSolver : public LinearSolver<double, SparseMatrixPolicy>
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP solve_function_resource_tracker_;
    void (*solve_function_)(const double*, double*, const double*, const double*);

   public:
    /// @brief Create a JITed linear solver for a given sparse matrix structure
    /// @param compiler JIT compiler
    /// @param matrix Block-diagonal sparse matrix to create solver for
    /// @param initial_value Initial value for decomposed triangular matrix elements
    JitLinearSolver(
        std::shared_ptr<JitCompiler> compiler,
        const SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
        double initial_value);

    ~JitLinearSolver();

    /// @brief Decompose the matrix into upper and lower triangular matrices and general JIT functions
    void Factor(SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix);

    /// @brief Solve for x in Ax = b
    template<template<class> class MatrixPolicy>
    void Solve(const MatrixPolicy<double>& b, MatrixPolicy<double>& x);

   private:
    /// @brief Generates the JIT-ed Solve function
    void GenerateSolveFunction();
  };

}  // namespace micm

#include "jit_linear_solver.inl"