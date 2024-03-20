// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <functional>
#include <micm/solver/cuda_linear_solver.cuh>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  template<typename T, template<class> class SparseMatrixPolicy, class LuDecompositionPolicy = CudaLuDecomposition>
  class CudaLinearSolver : public LinearSolver<T, SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    /// This is an instance of struct "LinearSolverParam" that holds
    ///   the constant data of "CudaLinearSolver" class on the device
    LinearSolverParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLinearSolver(){};

    /// This constructor takes two arguments: a sparse matrix and its values
    /// The base class here takes three arguments: the third argument is
    ///   a lamda function that creates an instance of LuDecompositionPolicy;
    ///   in this case, we will use the CudaLuDecomposition specified at line 15;
    ///   See line 62 of "linear_solver.inl" for more details about how
    ///   this lamda function works;
    CudaLinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value)
        : CudaLinearSolver<T, SparseMatrixPolicy, LuDecompositionPolicy>(
              matrix,
              initial_value,
              [&](const SparseMatrixPolicy<double>& m) -> LuDecompositionPolicy { return LuDecompositionPolicy(m); }){};

    CudaLinearSolver(
        const SparseMatrixPolicy<T>& matrix,
        T initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy<T>&)> create_lu_decomp)
        : LinearSolver<T, SparseMatrixPolicy, LuDecompositionPolicy>(matrix, initial_value, create_lu_decomp)
    {
      LinearSolverParam hoststruct;

      hoststruct.nLij_Lii_ = this->nLij_Lii_.data();
      hoststruct.Lij_yj_ = this->Lij_yj_.data();
      hoststruct.nUij_Uii_ = this->nUij_Uii_.data();
      hoststruct.Uij_xj_ = this->Uij_xj_.data();

      hoststruct.nLij_Lii_size_ = this->nLij_Lii_.size();
      hoststruct.Lij_yj_size_ = this->Lij_yj_.size();
      hoststruct.nUij_Uii_size_ = this->nUij_Uii_.size();
      hoststruct.Uij_xj_size_ = this->Uij_xj_.size();

      /// Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    }

    /// This is the destructor that will free the device memory of
    ///   the constant data from the class "CudaLinearSolver"
    ~CudaLinearSolver()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    template<template<class> class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy<T>> || VectorizableSparse<SparseMatrixPolicy<T>>) std::chrono::nanoseconds Solve(
        const MatrixPolicy<T>& b,
        MatrixPolicy<T>& x,
        SparseMatrixPolicy<T>& lower_matrix,
        SparseMatrixPolicy<T>& upper_matrix)
    {
      CudaSparseMatrixParam sparseMatrix;
      CudaMatrixParam denseMatrix;

      sparseMatrix.lower_matrix_ = lower_matrix.AsVector().data();
      sparseMatrix.lower_matrix_size_ = lower_matrix.AsVector().size();
      sparseMatrix.upper_matrix_ = upper_matrix.AsVector().data();
      sparseMatrix.upper_matrix_size_ = upper_matrix.AsVector().size();
      denseMatrix.b_ = b.AsVector().data();
      denseMatrix.x_ = x.AsVector().data();
      denseMatrix.b_size_ = b.AsVector().size();
      denseMatrix.x_size_ = x.AsVector().size();
      denseMatrix.n_grids_ = b.size();  // number of grids
      denseMatrix.b_column_counts_ = b[0].size();
      denseMatrix.x_column_counts_ = x[0].size();

      /// Call the "SolveKernelDriver" function that invokes the
      ///   CUDA kernel to perform the "solve" function on the device
      return micm::cuda::SolveKernelDriver(sparseMatrix, denseMatrix, this->devstruct_);
    };
  };
}  // namespace micm
