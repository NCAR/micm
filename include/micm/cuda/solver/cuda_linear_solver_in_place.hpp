// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_linear_solver_in_place.cuh>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/solver/linear_solver_in_place.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>

namespace micm
{
  template<class SparseMatrixPolicy, class LuDecompositionPolicy = CudaLuDecompositionMozartInPlace>
  class CudaLinearSolverInPlace : public LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    /// This is an instance of struct "LinearSolverInPlaceParam" that holds
    ///   the constant data of "CudaLinearSolverInPlace" class on the device
    LinearSolverInPlaceParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLinearSolverInPlace(){};

    CudaLinearSolverInPlace(const CudaLinearSolverInPlace&) = delete;
    CudaLinearSolverInPlace& operator=(const CudaLinearSolverInPlace&) = delete;
    CudaLinearSolverInPlace(CudaLinearSolverInPlace&& other)
        : LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };

    CudaLinearSolverInPlace& operator=(CudaLinearSolverInPlace&& other)
    {
      LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// This constructor takes two arguments: a sparse matrix and its values
    /// The base class here takes three arguments: the third argument is
    ///   a lamda function that creates an instance of LuDecompositionPolicy;
    ///   in this case, we will use the CudaLuDecompositionInPlace specified at line 13;
    ///   See line 17 of "linear_solver_in_place.inl" for more details about how
    ///   this lamda function works;
    CudaLinearSolverInPlace(const SparseMatrixPolicy& matrix, typename SparseMatrixPolicy::value_type initial_value)
        : CudaLinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>(
              matrix,
              initial_value,
              [&](const SparseMatrixPolicy& m) -> LuDecompositionPolicy { return LuDecompositionPolicy(m); }){};

    CudaLinearSolverInPlace(
        const SparseMatrixPolicy& matrix,
        typename SparseMatrixPolicy::value_type initial_value,
        const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)> create_lu_decomp)
        : LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>(matrix, initial_value, create_lu_decomp)
    {
      LinearSolverInPlaceParam hoststruct;

      hoststruct.nLij_ = this->nLij_.data();
      hoststruct.Lij_yj_ = this->Lij_yj_.data();
      hoststruct.nUij_Uii_ = this->nUij_Uii_.data();
      hoststruct.Uij_xj_ = this->Uij_xj_.data();

      hoststruct.nLij_size_ = this->nLij_.size();
      hoststruct.Lij_yj_size_ = this->Lij_yj_.size();
      hoststruct.nUij_Uii_size_ = this->nUij_Uii_.size();
      hoststruct.Uij_xj_size_ = this->Uij_xj_.size();

      /// Create the ALU matrix with all the fill-ins for the non-zero values
      auto ALU = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(matrix, 0, true);
      hoststruct.number_of_non_zeros_ = ALU.GroupSize() / SparseMatrixPolicy::GroupVectorSize();

      /// Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    }

    /// This is the destructor that will free the device memory of
    ///   the constant data from the class "CudaLinearSolverInPlace"
    ~CudaLinearSolverInPlace()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    template<class MatrixPolicy>
      requires(
          CudaMatrix<SparseMatrixPolicy> && CudaMatrix<MatrixPolicy> && VectorizableDense<MatrixPolicy> &&
          VectorizableSparse<SparseMatrixPolicy>)
    void Solve(MatrixPolicy& x, const SparseMatrixPolicy& ALU) const
    {
      auto x_param = x.AsDeviceParam();  // we need to update x so it can't be constant and must be an lvalue
      micm::cuda::SolveKernelDriver(x_param, ALU.AsDeviceParam(), this->devstruct_);
    };
  };
}  // namespace micm
