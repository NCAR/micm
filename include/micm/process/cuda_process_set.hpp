/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
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
    /// This is an instance of struct "ProcessSetParam" that holds
    ///   the constant data of "ProcessSet" class on the device
    ProcessSetParam devstruct_;

    CudaProcessSet() = default;
    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    CudaProcessSet(const std::vector<Process>& processes, const std::map<std::string, std::size_t>& variable_map);

    /// @brief Set the indexes for the elements of Jacobian matrix before we could copy it to the device;
    /// @brief this will override the "SetJacobianFlatIds" function from the "ProcessSet" class
    /// @param OrderingPolicy
    /// @param matrix
    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);

    template<template<class> typename MatrixPolicy>
    requires VectorizableDense<MatrixPolicy<double>>
    void AddForcingTerms(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& state_variables,
        MatrixPolicy<double>& forcing) const;

    template<class MatrixPolicy, class SparseMatrixPolicy>
    requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
    void SubtractJacobianTerms(
        const MatrixPolicy<double>& rate_constants,
        const MatrixPolicy<double>& state_variables,
        SparseMatrixPolicy<double>& jacobian) const;
  };

  inline CudaProcessSet::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::map<std::string, std::size_t>& variable_map)
      : ProcessSet(processes, variable_map)
  {
    /// Passing the class itself as an argument is not support by CUDA;
    /// Thus we generate a host struct first to save the pointers to
    ///   the actual data and size of each constant data member;

    /// Allocate host memory space for an object of type "ProcessSetParam"
    ProcessSetParam hoststruct;

    hoststruct.number_of_reactants_ = this->number_of_reactants_.data();
    hoststruct.reactant_ids_ = this->reactant_ids_.data();
    hoststruct.number_of_products_ = this->number_of_products_.data();
    hoststruct.product_ids_ = this->product_ids_.data();
    hoststruct.yields_ = this->yields_.data();
    hoststruct.jacobian_flat_ids_ = nullptr;

    hoststruct.number_of_reactants_size_ = this->number_of_reactants_.size();
    hoststruct.reactant_ids_size_ = this->reactant_ids_.size();
    hoststruct.number_of_products_size_ = this->number_of_products_.size();
    hoststruct.product_ids_size_ = this->product_ids_.size();
    hoststruct.yields_size_ = this->yields_.size();

    // Copy the data from host struct to device struct
    this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
  }

  template<typename OrderingPolicy>
  inline void CudaProcessSet::SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix)
  {
    /// This function sets the "jacobian_flat_ids_" member after the structure of Jacobian matrix is known
    ProcessSet::SetJacobianFlatIds(matrix);

    ProcessSetParam hoststruct;
    hoststruct.jacobian_flat_ids_ = this->jacobian_flat_ids_.data();
    hoststruct.jacobian_flat_ids_size_ = this->jacobian_flat_ids_.size();

    // Copy the data from host struct to device struct
    micm::cuda::CopyJacobiFlatId(hoststruct, this->devstruct_);
  }

  template<template<class> class MatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>>
  inline void CudaProcessSet::AddForcingTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      MatrixPolicy<double>& forcing) const
  {
    auto forcing_param = forcing.AsDeviceParam();  // we need to update forcing so it can't be constant and must be an lvalue
    micm::cuda::AddForcingTermsKernelDriver(
        rate_constants.AsDeviceParam(), state_variables.AsDeviceParam(), forcing_param, this->devstruct_);
  }

  template<class MatrixPolicy, class SparseMatrixPolicy>
  requires VectorizableDense<MatrixPolicy<double>> && VectorizableSparse<SparseMatrixPolicy<double>>
  inline void CudaProcessSet::SubtractJacobianTerms(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& state_variables,
      SparseMatrixPolicy<double>& jacobian) const
  {
    auto jacobian_param =
        jacobian.AsDeviceParam();  // we need to update jacobian so it can't be constant and must be an lvalue
    micm::cuda::SubtractJacobianTermsKernelDriver(
        rate_constants.AsDeviceParam(), state_variables.AsDeviceParam(), jacobian_param, this->devstruct_);
  }
}  // namespace micm
