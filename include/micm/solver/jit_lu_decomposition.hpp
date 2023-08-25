// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix optimized with JIT-compilation
  ///
  /// See LuDecomposition class description for algorithm details
  template<std::size_t L>
  class JitLuDecomposition : public LuDecomposition
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP decompose_function_resource_tracker_;
    void (*decompose_function_)(const double *, double *, double *);

   public:
    /// @brief Create a JITed LU decomposer for a given sparse matrix structure
    /// @param compiler JIT compiler
    /// @param matrix Sparse matrix to create LU decomposer for
    JitLuDecomposition(
        std::shared_ptr<JitCompiler> compiler,
        const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix);

    ~JitLuDecomposition();
  };

  template<std::size_t L>
  inline JitLuDecomposition<L>::JitLuDecomposition(
      std::shared_ptr<JitCompiler> compiler,
      const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix)
      : LuDecomposition(matrix),
        compiler_(compiler)
  {
    decompose_function_ = NULL;
  }

  template<std::size_t L>
  JitLuDecomposition<L>::~JitLuDecomposition()
  {
    if (decompose_function_ != NULL)
    {
      llvm::ExitOnError exit_on_error;
      exit_on_error(decompose_function_resource_tracker_->remove());
    }
  }

}  // namespace micm