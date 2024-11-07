// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <cublas_v2.h>

namespace micm
{
  /// @brief Construct a state variable for CUDA tests
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  struct CudaState : public State<DenseMatrixPolicy, SparseMatrixPolicy>
  {
   public:
    CudaState(const CudaState&) = delete;
    CudaState& operator=(const CudaState&) = delete;
    CudaState(CudaState&&) = default;
    CudaState& operator=(CudaState&&) = default;

    double* cuda_relative_tolerance_{ nullptr };

    ~CudaState()
    {
      CHECK_CUDA_ERROR(cudaFree(cuda_relative_tolerance_), "cudaFree");
    }

    /// @brief Constructor which takes the state dimension information as input
    /// @param parameters State dimension information
    CudaState(const StateParameters& parameters)
        : State<DenseMatrixPolicy, SparseMatrixPolicy>(parameters){};

    /// @brief Copy input variables to the device
    void SyncInputsToDevice() requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToDevice();
      this->rate_constants_.CopyToDevice();
      this->absolute_tolerance_.CopyToDevice();

      size_t relative_tolerance_bytes = sizeof(double);

      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &(cuda_relative_tolerance_), relative_tolerance_bytes, micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)),
          "cudaMalloc");
    }

    /// @brief Copy output variables to the host
    void SyncOutputsToHost() requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToHost();
    }
  };
}  // namespace micm