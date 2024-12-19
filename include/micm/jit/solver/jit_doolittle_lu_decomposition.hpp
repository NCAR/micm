// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/jit/jit_function.hpp>
#include <micm/solver/doolittle_lu_decomposition.hpp>
#include <micm/util/random_string.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief LU decomposer for SparseMatrix with vector-ordering optimized with JIT-compilation
  ///
  /// See DoolittleLuDecomposition class description for algorithm details
  /// The template parameter is the number of blocks (i.e. grid cells) in the block-diagonal matrix
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class JitDoolittleLuDecomposition : public DoolittleLuDecomposition
  {
    llvm::orc::ResourceTrackerSP decompose_function_resource_tracker_;
    using FuncPtr = void (*)(const double *, double *, double *);
    FuncPtr decompose_function_ = nullptr;

   public:
    JitDoolittleLuDecomposition(){};

    JitDoolittleLuDecomposition(const JitDoolittleLuDecomposition &) = delete;
    JitDoolittleLuDecomposition &operator=(const JitDoolittleLuDecomposition &) = delete;
    JitDoolittleLuDecomposition(JitDoolittleLuDecomposition &&other);
    JitDoolittleLuDecomposition &operator=(JitDoolittleLuDecomposition &&other);

    /// @brief Create a JITed LU decomposer for a given sparse matrix structure
    /// @param matrix Sparse matrix to create LU decomposer for
    JitDoolittleLuDecomposition(const SparseMatrix<double, SparseMatrixVectorOrdering<L>> &matrix);

    ~JitDoolittleLuDecomposition();

    /// @brief Create sparse L and U matrices for a given A matrix
    /// @param A Sparse matrix that will be decomposed
    /// @param lower The lower triangular matrix created by decomposition
    /// @param upper The upper triangular matrix created by decomposition
    template<class SparseMatrixPolicy>
    void Decompose(const SparseMatrixPolicy &A, SparseMatrixPolicy &lower, SparseMatrixPolicy &upper) const;

   private:
    /// @brief Generates a function to perform the LU decomposition for a specific matrix sparsity structure
    void GenerateDecomposeFunction();
  };

}  // namespace micm

#include "jit_doolittle_lu_decomposition.inl"