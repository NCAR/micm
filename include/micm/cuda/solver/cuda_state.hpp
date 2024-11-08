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

    CudaMatrixParam absolute_tolerance_param_;

    ~CudaState()
    {
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(absolute_tolerance_param_), "cudaFree");
      absolute_tolerance_param_.d_data_ = nullptr;
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

      absolute_tolerance_param_.number_of_elements_ = this->absolute_tolerance_.size();
      absolute_tolerance_param_.number_of_grid_cells_ = 1;

      CHECK_CUDA_ERROR(micm::cuda::MallocVector<double>(absolute_tolerance_param_, absolute_tolerance_param_.number_of_elements_), "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice<double>(absolute_tolerance_param_, this->absolute_tolerance_), "cudaMemcpyHostToDevice");
    }

    /// @brief Copy output variables to the host
    void SyncOutputsToHost() requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToHost();
    }
  };
}  // namespace micm