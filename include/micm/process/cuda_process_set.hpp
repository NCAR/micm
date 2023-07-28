
// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/process_set.hpp>

#ifdef USE_CUDA
#  include <micm/process/cuda_process_set.cuh>
#endif

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

#ifdef USE_CUDA
    template<template<class> typename MatrixPolicy>
    requires VectorizableDense<MatrixPolicy<double>>
    void AddForcingTerms(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& state_variables,
        MatrixPolicy<double>& forcing) const;
#endif
  };

  template<template<class> class MatrixPolicy>
  inline CudaProcessSet::CudaProcessSet(const std::vector<Process>& processes, const State<MatrixPolicy>& state)
      : ProcessSet(processes, state)
  {
  }

#ifdef USE_CUDA
  template<template<class> class MatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>>
  inline void CudaProcessSet::AddForcingTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      MatrixPolicy<double>& forcing) const
  {
    micm::cuda::AddForcingTerms_kernelSetup(
        rate_constants.AsVector().data(),
        state_variables.AsVector().data(),
        forcing.AsVector().data(),
        rate_constants.size(),
        rate_constants[0].size(),
        state_variables[0].size(),
        number_of_reactants_.data(),
        number_of_reactants_.size(),
        reactant_ids_.data(),
        reactant_ids_.size(),
        number_of_products_.data(),
        number_of_products_.size(),
        product_ids_.data(),
        product_ids_.size(),
        yields_.data(),
        yields_.size());
  }
#endif
}  // namespace micm