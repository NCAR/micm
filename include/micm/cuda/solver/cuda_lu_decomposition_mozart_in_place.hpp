// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>

namespace micm
{
  /// This CudaLuDecompositionMozartInPlace class inherits everything from the base class "LuDecompositionMozartInPlace"
  class CudaLuDecompositionMozartInPlace : public LuDecompositionMozartInPlace
  {
   public:
    /// This is an instance of struct "LuDecomposeMozartInPlaceParam" that holds
    ///   the constant data of "CudaLuDecompositionMozartInPlace" class on the device
    LuDecomposeMozartInPlaceParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLuDecompositionMozartInPlace(){};

    CudaLuDecompositionMozartInPlace(const CudaLuDecompositionMozartInPlace&) = delete;
    CudaLuDecompositionMozartInPlace& operator=(const CudaLuDecompositionMozartInPlace&) = delete;
    CudaLuDecompositionMozartInPlace(CudaLuDecompositionMozartInPlace&& other)
        : LuDecompositionMozartInPlace(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };

    CudaLuDecompositionMozartInPlace& operator=(CudaLuDecompositionMozartInPlace&& other)
    {
      LuDecompositionMozartInPlace::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// This is the overloaded constructor that takes one argument called "matrix";
    /// We need to specify the type (e.g., double, int, etc) and
    ///   ordering (e.g., vector-stored, non-vector-stored, etc) of the "matrix";
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    CudaLuDecompositionMozartInPlace(const SparseMatrixPolicy& matrix)
    {
      Initialize<SparseMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());

      /// Passing the class itself as an argument is not support by CUDA;
      /// Thus we generate a host struct first to save the pointers to
      ///   the actual data and size of each constant data member;

      /// Allocate host memory space for an object of type "LuDecomposeMozartInPlaceParam"
      LuDecomposeMozartInPlaceParam hoststruct;
      hoststruct.aii_nji_nki_ = this->aii_nji_nki_.data();
      hoststruct.aji_ = this->aji_.data();
      hoststruct.aik_njk_ = this->aik_njk_.data();
      hoststruct.ajk_aji_ = this->ajk_aji_.data();
      hoststruct.aii_nji_nki_size_ = this->aii_nji_nki_.size();
      hoststruct.aji_size_ = this->aji_.size();
      hoststruct.aik_njk_size_ = this->aik_njk_.size();
      hoststruct.ajk_aji_size_ = this->ajk_aji_.size();

      /// Create the ALU matrix with all the fill-ins for the non-zero values
      auto ALU = GetLUMatrix<SparseMatrixPolicy>(matrix, 0, true);
      hoststruct.number_of_non_zeros_ = ALU.GroupSize() / SparseMatrixPolicy::GroupVectorSize();

      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    /// This is destructor that will free the device memory of
    ///   the constant data from the class "CudaLuDecompositionMozartInPlace"
    ~CudaLuDecompositionMozartInPlace()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static CudaLuDecompositionMozartInPlace Create(const SparseMatrixPolicy& matrix)
    {
      CudaLuDecompositionMozartInPlace lu_decomp(matrix);
      return lu_decomp;
    }

    /// @brief This is the function to perform an LU decomposition on a given A matrix on the GPU
    /// @param ALU Sparse matrix to decompose (will be overwritten with L and U matrices)
    template<class SparseMatrixPolicy>
      requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(SparseMatrixPolicy& ALU) const;
  };

  template<class SparseMatrixPolicy>
    requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  void CudaLuDecompositionMozartInPlace::Decompose(SparseMatrixPolicy& ALU) const
  {
    auto ALU_param = ALU.AsDeviceParam();
    micm::cuda::DecomposeKernelDriver(ALU_param, this->devstruct_);
  }
}  // end of namespace micm
