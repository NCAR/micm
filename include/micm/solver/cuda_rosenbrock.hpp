// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process.hpp>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_rosenbrock.cuh>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/cuda_dense_matrix.hpp>
#include <micm/util/cuda_param.hpp>
#include <micm/util/cuda_sparse_matrix.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace micm
{

  template<class RatesPolicy, class LinearSolverPolicy>
  class CudaRosenbrockSolver : public AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>
  {
    ///@brief Default constructor
   public:
    /// This is an instance of struct "CudaRosenbrockSolverParam" that allocates
    ///   device memory of temporary variables and copy constant data member to device
    CudaRosenbrockSolverParam devstruct_;

    /// @brief Solver parameters typename
    using ParametersType = RosenbrockSolverParameters;

    CudaRosenbrockSolver(const CudaRosenbrockSolver&) = delete;
    CudaRosenbrockSolver& operator=(const CudaRosenbrockSolver&) = delete;
    CudaRosenbrockSolver(CudaRosenbrockSolver&& other)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(std::move(other)),
          devstruct_(std::move(other.devstruct_))
    {
      other.devstruct_.errors_input_ = nullptr;
      other.devstruct_.errors_output_ = nullptr;
      other.devstruct_.absolute_tolerance_ = nullptr;
      other.devstruct_.jacobian_diagonal_elements_ = nullptr;
    };

    CudaRosenbrockSolver& operator=(CudaRosenbrockSolver&& other)
    {
      RosenbrockSolver<RatesPolicy, LinearSolverPolicy>::operator=(std::move(other));
      devstruct_ = std::move(other.devstruct_);
      other.devstruct_.errors_input_ = nullptr;
      other.devstruct_.errors_output_ = nullptr;
      other.devstruct_.absolute_tolerance_ = nullptr;
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
        RosenbrockSolverParameters parameters,
        LinearSolverPolicy&& linear_solver,
        RatesPolicy&& rates,
        auto& jacobian)
        : AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>(parameters, std::move(linear_solver), std::move(rates), jacobian)
    {
      CudaRosenbrockSolverParam hoststruct;
      // jacobian.GroupVectorSize() is the same as the number of grid cells for the CUDA implementation
      // the absolute tolerance size is the same as the number of solved variables in one grid cell
      hoststruct.errors_size_ = jacobian.GroupVectorSize() * this->parameters_.absolute_tolerance_.size();
      hoststruct.jacobian_diagonal_elements_ = this->jacobian_diagonal_elements_.data();
      hoststruct.jacobian_diagonal_elements_size_ = this->jacobian_diagonal_elements_.size();
      hoststruct.absolute_tolerance_ = this->parameters_.absolute_tolerance_.data();
      hoststruct.absolute_tolerance_size_ = this->parameters_.absolute_tolerance_.size();
      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    /// This is the destructor that will free the device memory of
    ///   the constant data from the class "CudaRosenbrockSolver"
    ~CudaRosenbrockSolver()
    {
      /// Free the device memory allocated by the members of "devstruct_"
      micm::cuda::FreeConstData(this->devstruct_);
    };

    /// @brief Advances the given step over the specified time step on the GPU
    /// @param time_step Time [s] to advance the state by
    /// @return A struct containing results and a status code
    SolverResult Solve(double time_step, auto& state) noexcept
    {
      SolverResult result{};
      result = micm::AbstractRosenbrockSolver<RatesPolicy, LinearSolverPolicy, CudaRosenbrockSolver<RatesPolicy, LinearSolverPolicy>>::Solve(time_step, state);
  //    state.SyncOutputsToHost();
      return result;
    }

    /// @brief Computes [alpha * I - jacobian] on the GPU
    /// @tparam SparseMatrixPolicy
    /// @param jacobian Jacobian matrix
    /// @param alpha
    template<class SparseMatrixPolicy>
    void AlphaMinusJacobian(SparseMatrixPolicy& jacobian, const double& alpha) const
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
    double NormalizedError(const DenseMatrixPolicy& y_old, const DenseMatrixPolicy& y_new, const DenseMatrixPolicy& errors)
        const requires(CudaMatrix<DenseMatrixPolicy>&& VectorizableDense<DenseMatrixPolicy>)
    {
      return micm::cuda::NormalizedErrorDriver(
          y_old.AsDeviceParam(), y_new.AsDeviceParam(), errors.AsDeviceParam(), this->parameters_, this->devstruct_);
    }
  };  // end CudaRosenbrockSolver
}  // namespace micm