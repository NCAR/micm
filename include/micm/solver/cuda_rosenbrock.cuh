// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <micm/util/cuda_param.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <vector>
#include "cublas_v2.h"

namespace micm{
    namespace cuda{
      /// This is the function that will copy the constant data
      ///   members of class "CudaRosenbrockSolverParam" to the device
      ///   and allocate device memory for temporary variables;
      CudaRosenbrockSolverParam CopyConstData(CudaRosenbrockSolverParam& hoststruct);
   
      /// This is the function that will delete the constant data
      ///   members of class "CudaRosenbrockSolverParam" on the device
      void FreeConstData(CudaRosenbrockSolverParam& devstruct);

      /// @brief Compute alpha - J[i] for each element i at the diagnoal of matrix J
      /// @param h_jacobian sparse matrix on the host (will be replaced by the CudaSparseMatrix class)
      /// @param num_elements number of elements in the h_jacobian (will be replaced by the CudaSparseMatrix class)
      /// @param alpha scalar variable
      /// @param devstruct device struct including the locations of diagonal elements of the Jacobian matrix
      /// @return
      void AlphaMinusJacobianDriver(double* h_jacobian,
                                    const size_t num_elements,
                                    const double alpha,
                                    const CudaRosenbrockSolverParam& devstruct);

      /// @brief Computes the scaled norm of the matrix errors on the GPU; assume all the data are GPU resident already
      /// @param y_old_param matrix on the device
      /// @param y_new_param matrix on the device
      /// @param errors_param pre-computed error matrix from Rosenbrock solver
      /// @param ros_param struct of Rosenbrock solver parameters
      /// @param handle cublas handle
      /// @return the scaled norm of the matrix errors
      double NormalizedErrorDriver(const CudaVectorMatrixParam& y_old_param,
                                   const CudaVectorMatrixParam& y_new_param,
                                   const CudaVectorMatrixParam& errors_param,
                                   const RosenbrockSolverParameters& ros_param,
                                   cublasHandle_t handle,
                                   CudaRosenbrockSolverParam devstruct);
    }//end cuda
}//end micm