// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_set.hpp>
#include <micm/util/cuda_matrix_param.hpp>

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
    micm::CUDAMatrixParam matrixParam(rate_constants.AsVector()); 
    std::chrono::nanoseconds kernel_duration = micm::cuda::AddForcingTermsKernelDriver(
        matrixParam,
        state_variables.AsVector().data(),
        forcing.AsVector().data(),
        rate_constants.size(),
        rate_constants[0].size(),
        state_variables[0].size(),
        number_of_reactants_.data(),
        reactant_ids_.data(),
        reactant_ids_.size(),
        number_of_products_.data(),
        product_ids_.data(),
        product_ids_.size(),
        yields_.data(),
        yields_.size());
    return kernel_duration;  // time performance of kernel function
  }
  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
  inline std::chrono::nanoseconds CudaProcessSet::AddJacobianTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      SparseMatrixPolicy<double>& jacobian) const
  {
    std::chrono::nanoseconds kernel_duration = micm::cuda::AddJacobianTermsKernelDriver(
        rate_constants.AsVector().data(),
        state_variables.AsVector().data(),
        rate_constants.size(),      // n_grids
        rate_constants[0].size(),   // n_reactions
        state_variables[0].size(),  // n_species
        jacobian.AsVector().data(),
        jacobian.AsVector().size(),
        number_of_reactants_.data(),
        reactant_ids_.data(),
        reactant_ids_.size(),
        number_of_products_.data(),
        yields_.data(),
        yields_.size(),
        jacobian_flat_ids_.data(),
        jacobian_flat_ids_.size());
    return kernel_duration;  // time performance of kernel function
  }
}  // namespace micm
#endif