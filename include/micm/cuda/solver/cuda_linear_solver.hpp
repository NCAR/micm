// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_linear_solver.cuh>
#include <micm/cuda/solver/cuda_lu_decomposition.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>

#include <chrono>
#include <functional>

namespace micm
{
  template<class SparseMatrixPolicy, class LuDecompositionPolicy = CudaLuDecomposition>
  class CudaLinearSolver : public LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    /// This is an instance of struct "LinearSolverParam" that holds
    ///   the constant data of "CudaLinearSolver" class on the device
    LinearSolverParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLinearSolver(){};

    CudaLinearSolver(const CudaLinearSolver&) = delete;
    CudaLinearSolver& operator=(const CudaLinearSolver&) = delete;
    CudaLinearSolver(CudaLinearSolver&& other)
        : LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };

    CudaLinearSolver& operator=(CudaLinearSolver&& other)
    {
      LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// This constructor takes two arguments: a sparse matrix and its values
    /// The base class here takes three arguments: the third argument is
    ///   a lamda function that creates an instance of LuDecompositionPolicy;
    ///   in this case, we will use the CudaLuDecomposition specified at line 15;
    ///   See line 62 of "linear_solver.inl" for more details about how
    ///   this lamda function works;
    CudaLinearSolver(const SparseMatrixPolicy& matrix, typename SparseMatrixPolicy::value_type initial_value)
        : CudaLinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>(
              matrix,
              initial_value,
              [&](const SparseMatrixPolicy& m) -> LuDecompositionPolicy { return LuDecompositionPolicy(m); }){};

    CudaLinearSolver(
        const SparseMatrixPolicy& matrix,
        typename SparseMatrixPolicy::value_type initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)> create_lu_decomp)
        : LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>(matrix, initial_value, create_lu_decomp)
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

    template<class MatrixPolicy>
    requires(
        CudaMatrix<SparseMatrixPolicy>&& CudaMatrix<MatrixPolicy>&& VectorizableDense<MatrixPolicy>&& VectorizableSparse<
            SparseMatrixPolicy>) void Solve(MatrixPolicy& x, const SparseMatrixPolicy& L, const SparseMatrixPolicy& U) const
    {
      auto x_param = x.AsDeviceParam();  // we need to update x so it can't be constant and must be an lvalue
      micm::cuda::SolveKernelDriver(x_param, L.AsDeviceParam(), U.AsDeviceParam(), this->devstruct_);
    };
  };
}  // namespace micm
