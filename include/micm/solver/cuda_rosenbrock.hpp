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
    CudaRosenbrockSolver();

    CudaRosenbrockSolver(
        const System& system,
        const std::vector<Process>& processes,
        const RosenbrockSolverParameters& parameters)
        : RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(
              system,
              processes,
              parameters){};

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
              create_process_set){};

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
    /// @return
    double NormalizedError(MatrixPolicy<double>& y_old, 
                           MatrixPolicy<double>& y_new,
                           MatrixPolicy<double>& errors)
    {
      double* d_y_old = y_old.AsDeviceParam().d_data_;
      double* d_y_new = y_new.AsDeviceParam().d_data_;
      double* d_errors = errors.AsDeviceParam().d_data_;
      size_t num_elements = y_old.AsDeviceParam().num_elements_;
      double atol = this->parameters_.absolute_tolerance_;
      double rtol = this->parameters_.relative_tolerance_;

      if (y_old.AsDeviceParam().num_elements_ != y_new.AsDeviceParam().num_elements_ ||
          y_old.AsDeviceParam().num_elements_ != errors.AsDeviceParam().num_elements_)
      {
        throw std::runtime_error("The number of elements in y_old, y_new and errors must be the same.");
      }
      return micm::cuda::NormalizedErrorDriver(d_y_old, d_y_new, d_errors, num_elements, atol, rtol);
    }

  };  // end CudaRosenbrockSolver
}  // namespace micm