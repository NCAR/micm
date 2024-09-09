// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{
  /// @brief Construct a state variable for CUDA tests
  template<class TemporaryVariablesPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy>
  struct CudaState : public State<DenseMatrixPolicy, SparseMatrixPolicy>
  {
   public:
    CudaState(const CudaState&) = delete;
    CudaState& operator=(const CudaState&) = delete;
    CudaState(CudaState&&) = default;
    CudaState& operator=(CudaState&&) = default;

    /// @brief Constructor which takes the state dimension information as input
    /// @param parameters State dimension information
    CudaState(const StateParameters& parameters)
        : State<DenseMatrixPolicy, SparseMatrixPolicy>(parameters){};

    /// @brief Copy input variables to the device
    void SyncInputsToDevice() requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToDevice();
      this->rate_constants_.CopyToDevice();
    }

    /// @brief Copy output variables to the host
    void SyncOutputsToHost() requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToHost();
    }
  };
}  // namespace micm