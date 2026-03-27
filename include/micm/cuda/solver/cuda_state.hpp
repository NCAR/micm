// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{
  /// @brief Construct a state variable for CUDA tests
  template<
      class DenseMatrixPolicy = StandardDenseMatrix,
      class SparseMatrixPolicy = StandardSparseMatrix,
      class LuDecompositionPolicy = LuDecomposition>
  struct CudaState : public State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>
  {
   public:
    CudaState(const CudaState&) = delete;
    CudaState& operator=(const CudaState&) = delete;

    CudaMatrixParam absolute_tolerance_param_;

    CudaErrorParam errors_param_;
    CudaJacobianDiagonalElementsParam jacobian_diagonal_elements_param_;
    CudaConditionsParam conditions_param_;

    ~CudaState()
    {
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      CHECK_CUDA_ERROR(micm::cuda::FreeVector(absolute_tolerance_param_), "cudaFree");
      CHECK_CUDA_ERROR(micm::cuda::FreeArray(errors_param_.errors_input_), "cudaFree");
      CHECK_CUDA_ERROR(micm::cuda::FreeArray(errors_param_.errors_output_), "cudaFree");
      CHECK_CUDA_ERROR(micm::cuda::FreeArray(jacobian_diagonal_elements_param_.data_), "cudaFree");
      if (conditions_param_.d_temperature_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(conditions_param_.d_temperature_, stream), "cudaFree");
      if (conditions_param_.d_pressure_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(conditions_param_.d_pressure_, stream), "cudaFree");
      if (conditions_param_.d_air_density_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(conditions_param_.d_air_density_, stream), "cudaFree");
      if (conditions_param_.d_fixed_reactants_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(conditions_param_.d_fixed_reactants_, stream), "cudaFree");
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

      // If cuda vector length does not divide number of grid cells evenly,
      // we need to allocate GPU memory including the paddings for NormalizedError calculation.
      std::size_t cuda_rosenbrock_vector_length = DenseMatrixPolicy::GroupVectorSize();
      errors_param_.errors_size_ = parameters.number_of_species_ *
                                   ceil(static_cast<double>(number_of_grid_cells) / cuda_rosenbrock_vector_length) *
                                   cuda_rosenbrock_vector_length;

      auto diagonal_indices = this->jacobian_.DiagonalIndices(0);
      jacobian_diagonal_elements_param_.size_ = diagonal_indices.size();

      CHECK_CUDA_ERROR(
          micm::cuda::MallocVector<double>(absolute_tolerance_param_, absolute_tolerance_param_.number_of_elements_),
          "cudaMalloc");
      CHECK_CUDA_ERROR(
          micm::cuda::MallocArray<double>(errors_param_.errors_input_, errors_param_.errors_size_), "cudaMalloc");
      CHECK_CUDA_ERROR(
          micm::cuda::MallocArray<double>(errors_param_.errors_output_, errors_param_.errors_size_), "cudaMalloc");
      CHECK_CUDA_ERROR(
          micm::cuda::MallocArray<std::size_t>(
              jacobian_diagonal_elements_param_.data_, jacobian_diagonal_elements_param_.size_),
          "cudaMalloc");
      CHECK_CUDA_ERROR(micm::cuda::CopyToDevice<double>(absolute_tolerance_param_, atol), "cudaMemcpyHostToDevice");

      auto cuda_stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              jacobian_diagonal_elements_param_.data_,
              diagonal_indices.data(),
              sizeof(std::size_t) * diagonal_indices.size(),
              cudaMemcpyHostToDevice,
              cuda_stream),
          "cudaMemcpy");

      // Allocate persistent device buffers for conditions
      conditions_param_.number_of_grid_cells_ = number_of_grid_cells;
      conditions_param_.number_of_fixed_reactant_elements_ = this->rate_constants_.AsVector().size();
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&conditions_param_.d_temperature_, sizeof(double) * number_of_grid_cells, cuda_stream), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&conditions_param_.d_pressure_, sizeof(double) * number_of_grid_cells, cuda_stream), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&conditions_param_.d_air_density_, sizeof(double) * number_of_grid_cells, cuda_stream), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(
              &conditions_param_.d_fixed_reactants_,
              sizeof(double) * conditions_param_.number_of_fixed_reactant_elements_,
              cuda_stream),
          "cudaMalloc");
    };

    /// @brief Move constructor
    CudaState(CudaState&& other)
        : State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>(std::move(other))
    {
      absolute_tolerance_param_ = other.absolute_tolerance_param_;
      other.absolute_tolerance_param_.d_data_ = nullptr;
      errors_param_ = other.errors_param_;
      other.errors_param_ = {};
      jacobian_diagonal_elements_param_ = other.jacobian_diagonal_elements_param_;
      other.jacobian_diagonal_elements_param_ = {};
      conditions_param_ = other.conditions_param_;
      other.conditions_param_ = {};
    }

    /// @brief Move assignment operator
    CudaState& operator=(CudaState&& other)
    {
      if (this != &other)
      {
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::operator=(std::move(other));
        absolute_tolerance_param_ = other.absolute_tolerance_param_;
        other.absolute_tolerance_param_.d_data_ = nullptr;
        errors_param_ = other.errors_param_;
        other.errors_param_ = {};
        jacobian_diagonal_elements_param_ = other.jacobian_diagonal_elements_param_;
        other.jacobian_diagonal_elements_param_ = {};
        conditions_param_ = other.conditions_param_;
        other.conditions_param_ = {};
      }
      return *this;
    }

    void SetAbsoluteTolerances(const std::vector<double>& absoluteTolerance) override
    {
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>::SetAbsoluteTolerances(absoluteTolerance);
      CHECK_CUDA_ERROR(
          micm::cuda::CopyToDevice<double>(absolute_tolerance_param_, absoluteTolerance), "cudaMemcpyHostToDevice");
    }

    /// @brief Copy conditions (temperature, pressure, air_density) to the device
    void SyncConditionsToDevice()
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
    {
      const std::size_t number_of_grid_cells = conditions_param_.number_of_grid_cells_;
      std::vector<double> h_temperature(number_of_grid_cells);
      std::vector<double> h_pressure(number_of_grid_cells);
      std::vector<double> h_air_density(number_of_grid_cells);
      for (std::size_t i = 0; i < number_of_grid_cells; ++i)
      {
        h_temperature[i] = this->conditions_[i].temperature_;
        h_pressure[i] = this->conditions_[i].pressure_;
        h_air_density[i] = this->conditions_[i].air_density_;
      }
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              conditions_param_.d_temperature_, h_temperature.data(), sizeof(double) * number_of_grid_cells, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              conditions_param_.d_pressure_, h_pressure.data(), sizeof(double) * number_of_grid_cells, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              conditions_param_.d_air_density_, h_air_density.data(), sizeof(double) * number_of_grid_cells, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
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