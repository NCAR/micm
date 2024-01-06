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
  class CudaLuDecomposition : public LuDecomposition
  {
    public:
      /// This is an instance of struct "LuDecomposeConstDevice" that holds
      /// the constant data of "CudaLuDecomposition" class on the device
      LuDecomposeConstDevice* devptr = nullptr;

      /// This is the default constructor, taking no arguments;
      CudaLuDecomposition(){};

      /// This is the overloaded constructor that takes one argument called "matrix";
      /// We need to specify the type (e.g., double, int, etc) and 
      /// ordering (e.g., vector-stored, non-vector-stored, etc) of the "matrix";
      /// This CudaLuDecomposition class inherits everything from the base class "LuDecomposition";
      template<typename T, typename OrderingPolicy>
      CudaLuDecomposition(const SparseMatrix<T, OrderingPolicy>& matrix)
          : LuDecomposition(matrix)
      {
        /// Passing the class itself as an argument is not support by CUDA;
        /// Thus we generate a host struct first to save the pointers to 
        /// the value and size of each constant data member;

        /// Allocate memory space for ths object
        LuDecomposeConstHost* hostptr = new LuDecomposeConstHost;
        hostptr->niLU_         = this->niLU_.data();
        hostptr->do_aik_       = this->do_aik_.data();
        hostptr->aik_          = this->aik_.data();
        hostptr->uik_nkj_      = this->uik_nkj_.data();
        hostptr->lij_ujk_      = this->lij_ujk_.data();
        hostptr->do_aki_       = this->do_aki_.data();
        hostptr->aki_          = this->aki_.data();
        hostptr->lki_nkj_      = this->lki_nkj_.data();
        hostptr->lkj_uji_      = this->lkj_uji_.data();
        hostptr->uii_          = this->uii_.data();
        hostptr->niLU_size_    = this->niLU_.size();
        hostptr->do_aik_size_  = this->do_aik_.size();
        hostptr->aik_size_     = this->aik_.size();
        hostptr->uik_nkj_size_ = this->uik_nkj_.size();
        hostptr->lij_ujk_size_ = this->lij_ujk_.size();
        hostptr->do_aki_size_  = this->do_aki_.size();
        hostptr->aki_size_     = this->aki_.size();
        hostptr->lki_nkj_size_ = this->lki_nkj_.size();
        hostptr->lkj_uji_size_ = this->lkj_uji_.size();
        hostptr->uii_size_     = this->uii_.size(); 

        micm::cuda::CopyConstData(hostptr,this->devptr);

        /// Release the memory space for the hostptr
        delete hostptr->niLU_;
        delete hostptr->do_aik_;
        delete hostptr->aik_;
        delete hostptr->uik_nkj_;
        delete hostptr->lij_ujk_;
        delete hostptr->do_aki_;
        delete hostptr->aki_;
        delete hostptr->lki_nkj_;
        delete hostptr->lkj_uji_;
        delete hostptr->uii_;
        delete hostptr;
      };

      /// This is deconstructor that frees the device memory holding 
      /// the constant data from the CudaLuDecomposition class
      ~CudaLuDecomposition()
      {
        micm::cuda::FreeConstData(this->devptr);
      };

      /// This is the function to perform an LU decomposition on a given A matrix
      /// A is the sparse matrix to decompose
      /// L is the lower triangular matrix created by decomposition
      /// U is the upper triangular matrix created by decomposition
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
    CudaSparseMatrixParam sparseMatrix;
    sparseMatrix.A_ = A.AsVector().data();
    sparseMatrix.A_size_ = A.AsVector().size();
    sparseMatrix.L_ = L.AsVector().data();
    sparseMatrix.L_size_ = L.AsVector().size();
    sparseMatrix.U_ = U.AsVector().data();
    sparseMatrix.U_size_ = U.AsVector().size();
    sparseMatrix.n_grids_ = A.size();

    /// calling the DecomposeKernelDriver function that invokes the CUDA kernel to perform LU decomposition
    return micm::cuda::DecomposeKernelDriver(sparseMatrix, this->devptr);
  }
}  // end of namespace micm