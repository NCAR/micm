// Copyright (C) 2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/cuda_process_set.cuh>
#include <micm/process/process_set.hpp>
#include <micm/util/cuda_param.hpp>

namespace micm
{
  /// @brief A GPU-based implementation of ProcessSet
  class CudaProcessSet : public ProcessSet
  {
   public:
    CudaProcessSet() = default;
    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    CudaProcessSet(const std::vector<Process>& processes, const std::map<std::string, std::size_t>& variable_map);

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

  inline CudaProcessSet::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::map<std::string, std::size_t>& variable_map)
      : ProcessSet(processes, variable_map)
  {
  }

  template<template<class> class MatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>>
  inline std::chrono::nanoseconds CudaProcessSet::AddForcingTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      MatrixPolicy<double>& forcing) const
  {
    CudaMatrixParam matrix;
    matrix.rate_constants_ = rate_constants.AsVector().data();
    matrix.state_variables_ = state_variables.AsVector().data();
    matrix.forcing_ = forcing.AsVector().data();
    matrix.n_grids_ = rate_constants.size();
    matrix.n_reactions_ = rate_constants[0].size();
    matrix.n_species_ = state_variables[0].size();

    CudaProcessSetParam processSet;
    processSet.number_of_reactants_ = number_of_reactants_.data();
    processSet.reactant_ids_ = reactant_ids_.data();
    processSet.reactant_ids_size_ = reactant_ids_.size();
    processSet.number_of_products_ = number_of_products_.data();
    processSet.product_ids_ = product_ids_.data();
    processSet.product_ids_size_ = product_ids_.size();
    processSet.yields_ = yields_.data();
    processSet.yields_size_ = yields_.size();

    std::chrono::nanoseconds kernel_duration = micm::cuda::AddForcingTermsKernelDriver(matrix, processSet);
    return kernel_duration;  // time performance of kernel function
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
  inline std::chrono::nanoseconds CudaProcessSet::AddJacobianTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      SparseMatrixPolicy<double>& jacobian) const
  {
    CudaMatrixParam matrix;
    matrix.rate_constants_ = rate_constants.AsVector().data();
    matrix.state_variables_ = state_variables.AsVector().data();
    matrix.n_grids_ = rate_constants.size();
    matrix.n_reactions_ = rate_constants[0].size();
    matrix.n_species_ = state_variables[0].size();

    CudaSparseMatrixParam sparseMatrix;
    sparseMatrix.jacobian_ = jacobian.AsVector().data();
    sparseMatrix.jacobian_size_ = jacobian.AsVector().size();

    CudaProcessSetParam processSet;
    processSet.number_of_reactants_ = number_of_reactants_.data();
    processSet.reactant_ids_ = reactant_ids_.data();
    processSet.reactant_ids_size_ = reactant_ids_.size();
    processSet.number_of_products_ = number_of_products_.data();
    processSet.yields_ = yields_.data();
    processSet.yields_size_ = yields_.size();
    processSet.jacobian_flat_ids_ = jacobian_flat_ids_.data();
    processSet.jacobian_flat_ids_size_ = jacobian_flat_ids_.size();

    std::chrono::nanoseconds kernel_duration = micm::cuda::AddJacobianTermsKernelDriver(matrix, sparseMatrix, processSet);
    return kernel_duration;  // time performance of kernel function
  }
}  // namespace micm