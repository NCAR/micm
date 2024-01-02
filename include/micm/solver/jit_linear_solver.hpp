// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/jit_lu_decomposition.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/random_string.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief A general-use block-diagonal sparse-matrix linear solver with JIT-compiled optimizations
  ///
  /// See LinearSolver class description for algorithm details
  /// The template parameter is the number of blocks (i.e. grid cells) in the block-diagonal matrix
  template<std::size_t L, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy = JitLuDecomposition<L>>
  class JitLinearSolver : public LinearSolver<double, SparseMatrixPolicy, LuDecompositionPolicy>
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP solve_function_resource_tracker_;
    void (*solve_function_)(const double*, double*, const double*, const double*);

   public:
    JitLinearSolver(){};

    JitLinearSolver(const JitLinearSolver&) = delete;
    JitLinearSolver& operator=(const JitLinearSolver&) = delete;
    JitLinearSolver(JitLinearSolver&&);
    JitLinearSolver& operator=(JitLinearSolver&&);

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
    /// @param matrix Matrix that will be factored into lower and upper triangular matrices
    /// @param is_singular Flag that will be set to true if matrix is singular; false otherwise
    void Factor(
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& lower_matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& upper_matrix,
        bool& is_singular);

    /// @brief Decompose the matrix into upper and lower triangular matrices and general JIT functions
    /// @param matrix Matrix that will be factored into lower and upper triangular matrices
    void Factor(
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& lower_matrix,
        SparseMatrix<double, SparseMatrixVectorOrdering<L>>& upper_matrix);

    /// @brief Solve for x in Ax = b
    template<template<class> class MatrixPolicy>
    void Solve(
        const MatrixPolicy<double>& b,
        MatrixPolicy<double>& x,
        SparseMatrixPolicy<double>& lower_matrix,
        SparseMatrixPolicy<double>& upper_matrix);

   private:
    /// @brief Generates the JIT-ed Solve function
    void GenerateSolveFunction();
  };

}  // namespace micm

#include "jit_linear_solver.inl"