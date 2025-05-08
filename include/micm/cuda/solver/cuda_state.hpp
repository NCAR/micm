// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{
  /// @brief Construct a state variable for CUDA tests
  template<class DenseMatrixPolicy, class SparseMatrixPolicy, class LuDecompositionPolicy>
  struct CudaState : public State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    CudaState(const CudaState&) = delete;
    CudaState& operator=(const CudaState&) = delete;

    CudaMatrixParam absolute_tolerance_param_;

    ~CudaState()
    {
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(absolute_tolerance_param_), "cudaFree");
    }

    /// @brief Constructor which takes the state dimension information as input
    /// @param parameters State dimension information
    /// @param number_of_grid_cells Number of grid cells
    CudaState(const StateParameters& parameters, const std::size_t number_of_grid_cells)
        : State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>(parameters, number_of_grid_cells)
    {
      auto& atol = this->absolute_tolerance_;

      absolute_tolerance_param_.number_of_elements_ = atol.size();
      absolute_tolerance_param_.number_of_grid_cells_ = 1;

      CHECK_CUDA_ERROR(
          micm::cuda::MallocVector<double>(absolute_tolerance_param_, absolute_tolerance_param_.number_of_elements_),
          "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice<double>(absolute_tolerance_param_, atol), "cudaMemcpyHostToDevice");
    };

    /// @brief Move constructor
    CudaState(CudaState&& other)
        : State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>(std::move(other))
    {
      absolute_tolerance_param_ = other.absolute_tolerance_param_;
      other.absolute_tolerance_param_.d_data_ = nullptr;
    }

    /// @brief Move assignment operator
    CudaState& operator=(CudaState&& other)
    {
      if (this != &other)
      {
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::operator=(std::move(other));
        absolute_tolerance_param_ = other.absolute_tolerance_param_;
        other.absolute_tolerance_param_.d_data_ = nullptr;
      }
      return *this;
    }

    void SetAbsoluteTolerances(const std::vector<double>& absoluteTolerance) override
    {
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::SetAbsoluteTolerances(absoluteTolerance);
      CHECK_CUDA_ERROR(
          micm::cuda::CopyToDevice<double>(absolute_tolerance_param_, absoluteTolerance), "cudaMemcpyHostToDevice");
    }

    /// @brief Copy input variables to the device
    void SyncInputsToDevice()
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToDevice();
      this->rate_constants_.CopyToDevice();
    }

    /// @brief Copy output variables to the host
    void SyncOutputsToHost()
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
    {
      this->variables_.CopyToHost();
    }
  };
}  // namespace micm