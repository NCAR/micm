// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/random_string.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix with vector-ordering optimized with JIT-compilation
  ///
  /// See LuDecomposition class description for algorithm details
  /// The template parameter is the number of blocks (i.e. grid cells) in the block-diagonal matrix
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class JitLuDecomposition : public LuDecomposition
  {
    std::shared_ptr<JitCompiler> compiler_;
    llvm::orc::ResourceTrackerSP decompose_function_resource_tracker_;
    void (*decompose_function_)(const double *, double *, double *);

   public:
    JitLuDecomposition(){};

    JitLuDecomposition(const JitLuDecomposition &) = delete;
    JitLuDecomposition &operator=(const JitLuDecomposition &) = delete;
    JitLuDecomposition(JitLuDecomposition &&);
    JitLuDecomposition &operator=(JitLuDecomposition &&);

    /// @brief Create a JITed LU decomposer for a given sparse matrix structure
    /// @param compiler JIT compiler
    /// @param matrix Sparse matrix to create LU decomposer for
    JitLuDecomposition(
        std::shared_ptr<JitCompiler> compiler,
        const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix);

    ~JitLuDecomposition();

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @param lower The lower triangular matrix created by decomposition
    /// @param upper The upper triangular matrix created by decomposition
    /// @param is_singular Flag that will be set to true if A is singular; false otherwise
    template<typename T, template<class> class SparseMatrixPolicy>
    void Decompose(
        const SparseMatrixPolicy<T> &A,
        SparseMatrixPolicy<T> &lower,
        SparseMatrixPolicy<T> &upper,
        bool &is_singular) const;

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @param lower The lower triangular matrix created by decomposition
    /// @param upper The upper triangular matrix created by decomposition
    template<typename T, template<class> class SparseMatrixPolicy>
    void Decompose(const SparseMatrixPolicy<T> &A, SparseMatrixPolicy<T> &lower, SparseMatrixPolicy<T> &upper) const;

   private:
    /// @brief Generates a function to perform the LU decomposition for a specific matrix sparsity structure
    void GenerateDecomposeFunction();
  };

}  // namespace micm

#include "jit_lu_decomposition.inl"