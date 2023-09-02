// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_set.hpp>
#include <micm/util/cuda_param.hpp>

#ifdef USE_CUDA
#  include <micm/process/cuda_process_set.cuh>
#endif

#ifdef USE_CUDA
namespace micm
{
  /// @brief A GPU-based implementation of ProcessSet
  class CudaProcessSet : public ProcessSet
  {
   public:
    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param state Solver state
    template<template<class> class MatrixPolicy>
    CudaProcessSet(const std::vector<Process>& processes, const State<MatrixPolicy>& state);

    template<template<class> typename MatrixPolicy>
    requires VectorizableDense<MatrixPolicy<double>> std::chrono::nanoseconds AddForcingTerms(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& state_variables,
        MatrixPolicy<double>& forcing)
    const;

    template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
    requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
        std::chrono::nanoseconds AddJacobianTerms(
            const MatrixPolicy<double>& rate_constants,
            const MatrixPolicy<double>& state_variables,
            SparseMatrixPolicy<double>& jacobian)
    const;
  };

  template<template<class> class MatrixPolicy>
  inline CudaProcessSet::CudaProcessSet(const std::vector<Process>& processes, const State<MatrixPolicy>& state)
      : ProcessSet(processes, state)
  {
  }

  template<template<class> class MatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>>
  inline std::chrono::nanoseconds CudaProcessSet::AddForcingTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      MatrixPolicy<double>& forcing) const
  {
    CUDAMatrixParam matrixParam;
    matrixParam.rate_constants = rate_constants.AsVector().data(); 
    matrixParam.state_variables = state_variables.AsVector().data(); 
    matrixParam.forcing = forcing.AsVector().data(); 
    matrixParam.n_grids = rate_constants.size(); 
    matrixParam.n_reactions = rate_constants[0].size(); 
    matrixParam.n_species = state_variables[0].size();

    CUDAProcessSetParam processSet; 
    processSet.number_of_reactants = number_of_reactants_.data(); 
    processSet.reactant_ids = reactant_ids_.data(); 
    processSet.reactant_ids_size = reactant_ids_.size(); 
    processSet.number_of_products = number_of_products_.data(); 
    processSet.product_ids = product_ids_.data(); 
    processSet.product_ids_size = product_ids_.size(); 
    processSet.yields = yields_.data(); 
    processSet.yields_size = yields_.size(); 

    std::chrono::nanoseconds kernel_duration = micm::cuda::AddForcingTermsKernelDriver(
        matrixParam,
        processSet);
    return kernel_duration;  // time performance of kernel function
  }
  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
  inline std::chrono::nanoseconds CudaProcessSet::AddJacobianTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      SparseMatrixPolicy<double>& jacobian) const
  {
    CUDAMatrixParam matrixParam; 
    matrixParam.rate_constants = rate_constants.AsVector().data(); 
    matrixParam.state_variables = state_variables.AsVector().data(); 
    matrixParam.n_grids = rate_constants.size(); 
    matrixParam.n_reactions = rate_constants[0].size(); 
    matrixParam.n_species = state_variables[0].size(); 
    
    CUDASparseMatrixParam sparseMatrix; 
    sparseMatrix.jacobian = jacobian.AsVector().data(); 
    sparseMatrix.jacobian_size = jacobian.AsVector().size(); 
    
    CUDAProcessSetParam processSet; 
    processSet.number_of_reactants = number_of_reactants_.data(); 
    processSet.reactant_ids = reactant_ids_.data(); 
    processSet.reactant_ids_size = reactant_ids_.size(); 
    processSet.number_of_products = number_of_products_.data(); 
    processSet.yields = yields_.data(); 
    processSet.yields_size = yields_.size(); 
    processSet.jacobian_flat_ids = jacobian_flat_ids_.data(); 
    processSet.jacobian_flat_ids_size = jacobian_flat_ids_.size(); 
    
    std::chrono::nanoseconds kernel_duration = micm::cuda::AddJacobianTermsKernelDriver(
        matrixParam, 
        sparseMatrix, 
        processSet);
    return kernel_duration;  // time performance of kernel function
  }
}  // namespace micm
#endif