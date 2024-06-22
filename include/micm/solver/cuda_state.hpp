// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/cuda_dense_matrix.hpp>

namespace micm
{
  /// @brief Construct a state variable for CUDA tests
  template<class DenseMatrixPolicy = StandardDenseMatrix, class SparseMatrixPolicy = StandardSparseMatrix>
  struct CudaState : public State<DenseMatrixPolicy, SparseMatrixPolicy>
  {
    public:
      CudaState(const CudaState&) = delete;
      CudaState& operator=(const CudaState&) = delete;
      CudaState(CudaState&&) = default;
      CudaState& operator=(CudaState&&) = default;

      /// @brief Copy input variables to the device
      template<class DenseMatrixPolicy, class SparseMatrixPolicy>
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>) 
      void SyncInputsToDevice()
      {
         variables_.CopyToDevice();
         rate_constants_.CopyToDevice();
      }

      /// @brief Copy output variables to the host
      template<class DenseMatrixPolicy, class SparseMatrixPolicy>
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>) 
      void SyncOutputsToHost()
      {
        variables_.CopyToHost();
      }
  };
}  // namespace micm