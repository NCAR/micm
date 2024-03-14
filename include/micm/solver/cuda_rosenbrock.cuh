// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <chrono>
#include <micm/util/cuda_param.hpp>
#include <vector>
#include "cublas_v2.h"

namespace micm{
    namespace cuda{
      /// @brief Compute alpha - J[i] for each element i at the diagnoal of matrix J
      /// @param sparseMatrix struct of sparse matrix (will be replaced by the CudaSparseMatrix class)
      /// @param jacobian_diagonal_elements location id of the diagonal elements 
      /// @param alpha scalar variable
      /// @return
      std::chrono::nanoseconds AlphaMinusJacobianDriver(
                               CudaSparseMatrixParam& sparseMatrix,
                               const std::vector<size_t> jacobian_diagonal_elements,
                               double alpha);

      /// @brief Computes the scaled norm of the vector errors on the GPU; assume all the data are GPU resident already
      /// @param d_y_old device pointer to one vector
      /// @param d_y_new deviec pointer to another vector
      /// @param d_errors device pointer to the pre-computed error vector from Rosenbrock solver
      /// @param num_elements number of elements in the vectors
      /// @param atol absolute tolerance to be used in the error computation for a particular Rosenbrock solver
      /// @param rtol relative tolerance to be used in the error computation for a particular Rosenbrock solver
      /// @param handle cublas handle
      /// @return the scaled norm of the vector errors
      double NormalizedErrorDriver(double* d_y_old, double* d_y_new, 
                                   double* d_errors, const size_t num_elements,
                                   const double atol, const double rtol,
                                   cublasHandle_t handle);
    }//end cuda
}//end micm