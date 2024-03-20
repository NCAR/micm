// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <chrono>
#include <micm/solver/cuda_lu_decomposition.cuh>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/cuda_param.hpp>
#include <stdexcept>

namespace micm
{
  /// This CudaLuDecomposition class inherits everything from the base class "LuDecomposition"
  class CudaLuDecomposition : public LuDecomposition
  {
   public:
    /// This is an instance of struct "LuDecomposeParam" that holds
    ///   the constant data of "CudaLuDecomposition" class on the device
    LuDecomposeParam devstruct_;

    /// This is the default constructor, taking no arguments;
    CudaLuDecomposition(){};

    /// This is the overloaded constructor that takes one argument called "matrix";
    /// We need to specify the type (e.g., double, int, etc) and
    ///   ordering (e.g., vector-stored, non-vector-stored, etc) of the "matrix";
    template<typename T, typename OrderingPolicy>
    CudaLuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix)
        : LuDecomposition(matrix)
    {
      /// Passing the class itself as an argument is not support by CUDA;
      /// Thus we generate a host struct first to save the pointers to
      ///   the actual data and size of each constant data member;

      /// Allocate host memory space for an object of type "LuDecomposeParam"
      LuDecomposeParam hoststruct;
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
    ///   the constant data from the class "CudaLuDecomposition"
    ~CudaLuDecomposition()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// This is the function to perform an LU decomposition on a given A matrix
    ///   A is the sparse matrix to decompose
    ///   L is the lower triangular matrix created by decomposition
    ///   U is the upper triangular matrix created by decomposition
    template<typename T, template<class> typename SparseMatrixPolicy>
    requires VectorizableSparse<SparseMatrixPolicy<T>> std::chrono::nanoseconds
    Decompose(const SparseMatrixPolicy<T>& A, SparseMatrixPolicy<T>& L, SparseMatrixPolicy<T>& U)
    const;
  };

  template<typename T, template<class> class SparseMatrixPolicy>
  requires(VectorizableSparse<SparseMatrixPolicy<T>>) std::chrono::nanoseconds
      CudaLuDecomposition::Decompose(const SparseMatrixPolicy<T>& A, SparseMatrixPolicy<T>& L, SparseMatrixPolicy<T>& U)
  const
  {
    /// Once the CudaMatrix class is generated, we won't need the following lines any more;
    CudaSparseMatrixParam sparseMatrix;
    sparseMatrix.A_ = A.AsVector().data();
    sparseMatrix.A_size_ = A.AsVector().size();
    sparseMatrix.L_ = L.AsVector().data();
    sparseMatrix.L_size_ = L.AsVector().size();
    sparseMatrix.U_ = U.AsVector().data();
    sparseMatrix.U_size_ = U.AsVector().size();
    sparseMatrix.n_grids_ = A.size();

    /// Call the "DecomposeKernelDriver" function that invokes the
    ///   CUDA kernel to perform LU decomposition on the device
    return micm::cuda::DecomposeKernelDriver(sparseMatrix, this->devstruct_);
  }
}  // end of namespace micm
