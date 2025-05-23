// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/solver/cuda_linear_solver_in_place.hpp>
#include <micm/cuda/solver/cuda_rosenbrock.cuh>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

namespace micm
{
  struct CudaRosenbrockSolverParameters;

  template<class RatesPolicy, class LinearSolverPolicy>
  class CudaRosenbrockSolver : public AbstractRosenbrockSolver<
                                   RatesPolicy,
                                   LinearSolverPolicy,
                                   CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>
  {
    ///@brief Default constructor
   public:
    /// This is an instance of struct "CudaRosenbrockSolverParam" that allocates
    ///   device memory of temporary variables and copy constant data member to device
    CudaRosenbrockSolverParam devstruct_;

    /// @brief Solver parameters typename
    using ParametersType = CudaRosenbrockSolverParameters;

    CudaRosenbrockSolver(const CudaRosenbrockSolver&) = delete;
    CudaRosenbrockSolver& operator=(const CudaRosenbrockSolver&) = delete;
    CudaRosenbrockSolver(CudaRosenbrockSolver&& other)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(
              std::move(other)),
          devstruct_(std::move(other.devstruct_))
    {
      other.devstruct_.errors_input_ = nullptr;
      other.devstruct_.errors_output_ = nullptr;
      other.devstruct_.jacobian_diagonal_elements_ = nullptr;
    };

    CudaRosenbrockSolver& operator=(CudaRosenbrockSolver&& other)
    {
      RosenbrockSolver<RatesPolicy, LinearSolverPolicy>::operator=(std::move(other));
      devstruct_ = std::move(other.devstruct_);
      other.devstruct_.errors_input_ = nullptr;
      other.devstruct_.errors_output_ = nullptr;
      other.devstruct_.jacobian_diagonal_elements_ = nullptr;
      return *this;
    };

    /// @brief Default constructor
    CudaRosenbrockSolver()
    {
      devstruct_.errors_input_ = nullptr;
      devstruct_.errors_output_ = nullptr;
    };

    /// @brief Builds a CUDA Rosenbrock solver for the given system and solver parameters
    /// @param parameters Solver parameters
    /// @param linear_solver Linear solver
    /// @param rates Rates calculator
    /// @param jacobian Jacobian matrix
    CudaRosenbrockSolver(
        LinearSolverPolicy&& linear_solver,
        RatesPolicy&& rates,
        auto& jacobian,
        const size_t number_of_species)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(
              std::move(linear_solver),
              std::move(rates),
              jacobian,
              number_of_species)
    {
      CudaRosenbrockSolverParam hoststruct;
      // jacobian.GroupVectorSize() is the same as the number of grid cells for the CUDA implementation
      // the absolute tolerance size is the same as the number of solved variables in one grid cell
      hoststruct.errors_size_ = jacobian.GroupVectorSize() * number_of_species;
      auto diagonal_indices = jacobian.DiagonalIndices(0);
      hoststruct.jacobian_diagonal_elements_ = diagonal_indices.data();
      hoststruct.jacobian_diagonal_elements_size_ = diagonal_indices.size();
      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    /// This is the destructor that will free the device memory of
    ///   the constant data from the class "CudaRosenbrockSolver"
    ~CudaRosenbrockSolver()
    {
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// @brief Computes [alpha * I - jacobian] on the GPU
    /// @tparam SparseMatrixPolicy
    /// @param jacobian Jacobian matrix
    /// @param jacobian_diagonal_elements Diagonal elements of the Jacobian matrix, not used
    /// @param alpha
    template<class SparseMatrixPolicy>
    void AlphaMinusJacobian(
        SparseMatrixPolicy& jacobian,
        const std::vector<std::size_t>& jacobian_diagonal_elements,
        const double& alpha) const
      requires(CudaMatrix<SparseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
    {
      auto jacobian_param =
          jacobian.AsDeviceParam();  // we need to update jacobian so it can't be constant and must be an lvalue
      micm::cuda::AlphaMinusJacobianDriver(jacobian_param, alpha, this->devstruct_);
    }

    /// @brief Computes the scaled norm of the vector errors on the GPU; assume all the data are GPU resident already
    /// @tparam DenseMatrixPolicy
    /// @param y_old the original vector
    /// @param y_new the new vector
    /// @param errors The computed errors
    /// @return The scaled norm of the errors
    template<class DenseMatrixPolicy>
    double NormalizedError(
        const DenseMatrixPolicy& y_old,
        const DenseMatrixPolicy& y_new,
        const DenseMatrixPolicy& errors,
        auto& state) const
      requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
    {
      return micm::cuda::NormalizedErrorDriver(
          y_old.AsDeviceParam(),
          y_new.AsDeviceParam(),
          errors.AsDeviceParam(),
          state.absolute_tolerance_param_,
          state.relative_tolerance_,
          this->devstruct_);
    }
  };  // end CudaRosenbrockSolver
}  // namespace micm