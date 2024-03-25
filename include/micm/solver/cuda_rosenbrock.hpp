#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/cuda_process_set.hpp>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_rosenbrock.cuh>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/cuda_param.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <string>
#include <vector>

namespace micm
{

  template<
      template<class> class MatrixPolicy = Matrix,
      template<class> class SparseMatrixPolicy = StandardSparseMatrix,
      class LinearSolverPolicy = CudaLinearSolver<double, SparseMatrixPolicy>,
      class ProcessSetPolicy = CudaProcessSet>

  class CudaRosenbrockSolver
      : public RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
  {
    ///@brief Default constructor
   public:
    /// This is an instance of struct "CudaRosenbrockSolverParam" that allocates
    ///   device memory of temporary variables and copy constant data member to device
    CudaRosenbrockSolverParam devstruct_;

    CudaRosenbrockSolver()
    {
      devstruct_.errors_input_ = nullptr;
      devstruct_.errors_output_ = nullptr;
    };

    CudaRosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters)
        : RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(
              system,
              processes,
              parameters)
    {
      CudaRosenbrockSolverParam hoststruct;
      //      hoststruct.jacobian_diagonal_elements_ = state_parameters_.jacobian_diagonal_elements_;
      //      hoststruct.jacobian_diagonal_elements_size_ = state_parameters_.jacobian_diagonal_elements_.size();
      hoststruct.errors_size_ = parameters.number_of_grid_cells_ * system.StateSize();
      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    CudaRosenbrockSolver(
        const System& system,
        const std::vector<Process> processes,
        const RosenbrockSolverParameters& parameters,
        const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
        const std::function<ProcessSetPolicy(const std::vector<Process>&, const std::map<std::string, std::size_t>&)>
            create_process_set)
        : RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(
              system,
              processes,
              parameters,
              create_linear_solver,
              create_process_set)
    {
      CudaRosenbrockSolverParam hoststruct;
      //      hoststruct.jacobian_diagonal_elements_ = state_parameters_.jacobian_diagonal_elements_;
      //      hoststruct.jacobian_diagonal_elements_size_ = state_parameters_.jacobian_diagonal_elements_.size();
      hoststruct.errors_size_ = parameters.number_of_grid_cells_ * system.StateSize();
      // Copy the data from host struct to device struct
      this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
    };

    std::chrono::nanoseconds AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, double alpha) const requires
        VectorizableSparse<SparseMatrixPolicy<double>>
    {
      for (auto& element : jacobian.AsVector())
        element = -element;
      CudaSparseMatrixParam sparseMatrix;
      sparseMatrix.jacobian_ = jacobian.AsVector().data();
      sparseMatrix.jacobian_size_ = jacobian.AsVector().size();
      sparseMatrix.n_grids_ = jacobian.size();

      return micm::cuda::AlphaMinusJacobianDriver(sparseMatrix, this->state_parameters_.jacobian_diagonal_elements_, alpha);
    }

    /// @brief Computes the scaled norm of the vector errors on the GPU; assume all the data are GPU resident already
    /// @param y_old the original vector
    /// @param y_new the new vector
    /// @param errors The computed errors
    /// @return The scaled norm of the errors
    double NormalizedError(
        const MatrixPolicy<double>& y_old,
        const MatrixPolicy<double>& y_new,
        const MatrixPolicy<double>& errors) const
    {
      // At this point, it does not matter which handle we use; may revisit it when we have a multi-node-multi-GPU test
      return micm::cuda::NormalizedErrorDriver(
          y_old.AsDeviceParam(),
          y_new.AsDeviceParam(),
          errors.AsDeviceParam(),
          this->parameters_,
          errors.AsCublasHandle(),
          this->devstruct_);
    }
  };  // end CudaRosenbrockSolver
}  // namespace micm