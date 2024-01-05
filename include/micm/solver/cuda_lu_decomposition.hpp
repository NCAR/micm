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
        micm::cuda::CopyConstData(this,this->devptr);
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
