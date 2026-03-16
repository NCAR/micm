// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_linear_solver_in_place.cuh>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/solver/linear_solver_in_place.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>

#include <fstream>
#include <iostream>

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

      // Write const index arrays to a text file
      {
        std::ofstream outfile("linear_solver_in_place_const_arrays.txt");
        if (outfile.is_open())
        {
          outfile << "nLij_size: " << this->nLij_.size() << "\n";
          outfile << "nLij:";
          for (const auto& v : this->nLij_)
            outfile << " " << v;
          outfile << "\n";

          outfile << "Lij_yj_size: " << this->Lij_yj_.size() << "\n";
          outfile << "Lij_yj_first:";
          for (const auto& p : this->Lij_yj_)
            outfile << " " << p.first;
          outfile << "\n";
          outfile << "Lij_yj_second:";
          for (const auto& p : this->Lij_yj_)
            outfile << " " << p.second;
          outfile << "\n";

          outfile << "nUij_Uii_size: " << this->nUij_Uii_.size() << "\n";
          outfile << "nUij_Uii_first:";
          for (const auto& p : this->nUij_Uii_)
            outfile << " " << p.first;
          outfile << "\n";
          outfile << "nUij_Uii_second:";
          for (const auto& p : this->nUij_Uii_)
            outfile << " " << p.second;
          outfile << "\n";

          outfile << "Uij_xj_size: " << this->Uij_xj_.size() << "\n";
          outfile << "Uij_xj_first:";
          for (const auto& p : this->Uij_xj_)
            outfile << " " << p.first;
          outfile << "\n";
          outfile << "Uij_xj_second:";
          for (const auto& p : this->Uij_xj_)
            outfile << " " << p.second;
          outfile << "\n";

          outfile << "number_of_non_zeros: " << hoststruct.number_of_non_zeros_ << "\n";

          outfile.close();
          std::cout << "Wrote linear solver const arrays to: linear_solver_in_place_const_arrays.txt" << std::endl;
        }
      }
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
