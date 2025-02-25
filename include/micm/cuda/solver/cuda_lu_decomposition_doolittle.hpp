// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_lu_decomposition_doolittle.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition_doolittle.hpp>

namespace micm
{
  /// This CudaLuDecompositionDoolittle class inherits everything from the base class "LuDecomposition"
  class CudaLuDecompositionDoolittle : public LuDecompositionDoolittle
  {
   public:
    /// This is an instance of struct "LuDecomposeDoolittleParam" that holds
    ///   the constant data of "CudaLuDecompositionDoolittle" class on the device
    LuDecomposeDoolittleParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLuDecompositionDoolittle(){};

    CudaLuDecompositionDoolittle(const CudaLuDecompositionDoolittle&) = delete;
    CudaLuDecompositionDoolittle& operator=(const CudaLuDecompositionDoolittle&) = delete;
    CudaLuDecompositionDoolittle(CudaLuDecompositionDoolittle&& other)
        : LuDecompositionDoolittle(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };

    CudaLuDecompositionDoolittle& operator=(CudaLuDecompositionDoolittle&& other)
    {
      LuDecompositionDoolittle::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// This is the overloaded constructor that takes one argument called "matrix";
    /// We need to specify the type (e.g., double, int, etc) and
    ///   ordering (e.g., vector-stored, non-vector-stored, etc) of the "matrix";
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
    CudaLuDecompositionDoolittle(const SparseMatrixPolicy& matrix)
    {
      Initialize<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(matrix, typename SparseMatrixPolicy::value_type());

      /// Passing the class itself as an argument is not support by CUDA;
      /// Thus we generate a host struct first to save the pointers to
      ///   the actual data and size of each constant data member;

      /// Allocate host memory space for an object of type "LuDecomposeDoolittleParam"
      LuDecomposeDoolittleParam hoststruct;
      hoststruct.niLU_ = this->niLU_.data();
      hoststruct.do_aik_ = this->do_aik_.data();
      hoststruct.aik_ = this->aik_.data();
      hoststruct.uik_nkj_ = this->uik_nkj_.data();
      hoststruct.lij_ujk_ = this->lij_ujk_.data();
      hoststruct.do_aki_ = this->do_aki_.data();
      hoststruct.aki_ = this->aki_.data();
      hoststruct.lki_nkj_ = this->lki_nkj_.data();
      hoststruct.lkj_uji_ = this->lkj_uji_.data();
      hoststruct.uii_ = this->uii_.data();
      hoststruct.niLU_size_ = this->niLU_.size();
      hoststruct.do_aik_size_ = this->do_aik_.size();
      hoststruct.aik_size_ = this->aik_.size();
      hoststruct.uik_nkj_size_ = this->uik_nkj_.size();
      hoststruct.lij_ujk_size_ = this->lij_ujk_.size();
      hoststruct.do_aki_size_ = this->do_aki_.size();
      hoststruct.aki_size_ = this->aki_.size();
      hoststruct.lki_nkj_size_ = this->lki_nkj_.size();
      hoststruct.lkj_uji_size_ = this->lkj_uji_.size();
      hoststruct.uii_size_ = this->uii_.size();

      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    /// This is destructor that will free the device memory of
    ///   the constant data from the class "CudaLuDecompositionDoolittle"
    ~CudaLuDecompositionDoolittle()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// @brief Create an LU decomposition algorithm for a given sparse matrix policy
    /// @param matrix Sparse matrix
    template<class SparseMatrixPolicy, class LMatrixPolicy = SparseMatrixPolicy, class UMatrixPolicy = SparseMatrixPolicy>
      requires(SparseMatrixConcept<SparseMatrixPolicy>)
    static CudaLuDecompositionDoolittle Create(const SparseMatrixPolicy& matrix)
    {
      static_assert(
          std::is_same_v<SparseMatrixPolicy, LMatrixPolicy>,
          "SparseMatrixPolicy must be the same as LMatrixPolicy for CUDA LU decomposition");
      static_assert(
          std::is_same_v<SparseMatrixPolicy, UMatrixPolicy>,
          "SparseMatrixPolicy must be the same as UMatrixPolicy for CUDA LU decomposition");
      CudaLuDecompositionDoolittle lu_decomp(matrix);
      return lu_decomp;
    }

    /// @brief This is the function to perform an LU decomposition on a given A matrix on the GPU
    /// @param A is the sparse matrix to decompose
    /// @param L is the lower triangular matrix created by decomposition
    /// @param U is the upper triangular matrix created by decomposition
    template<class SparseMatrixPolicy>
      requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    void Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const;
  };

  template<class SparseMatrixPolicy>
    requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  void CudaLuDecompositionDoolittle::Decompose(const SparseMatrixPolicy& A, auto& L, auto& U) const
  {
    auto L_param = L.AsDeviceParam();  // we need to update lower matrix so it can't be constant and must be an lvalue
    auto U_param = U.AsDeviceParam();  // we need to update upper matrix so it can't be constant and must be an lvalue
    micm::cuda::DecomposeKernelDriver(A.AsDeviceParam(), L_param, U_param, this->devstruct_);
  }
}  // end of namespace micm
