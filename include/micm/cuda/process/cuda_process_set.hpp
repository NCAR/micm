// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/process/cuda_process_set.cuh>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/process_set.hpp>

namespace micm
{
  /// @brief A GPU-based implementation of ProcessSet
  /// @tparam DenseMatrixPolicy Policy for dense matrices (must satisfy CudaMatrix concept)
  /// @tparam SparseMatrixPolicy Policy for sparse matrices (must satisfy CudaMatrix concept)
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  class CudaProcessSet : public ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>
  {
   public:
    /// This is an instance of struct "ProcessSetParam" that holds
    ///   the constant data of "ProcessSet" class on the device
    ProcessSetParam devstruct_;

    CudaProcessSet() = default;

    CudaProcessSet(const CudaProcessSet&) = delete;
    CudaProcessSet& operator=(const CudaProcessSet&) = delete;
    CudaProcessSet(CudaProcessSet&& other)
        : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(std::move(other))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };
    CudaProcessSet& operator=(CudaProcessSet&& other)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      return *this;
    };

    /// @brief Create a process set calculator for a given set of processes
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    CudaProcessSet(const std::vector<Process>& processes, const std::unordered_map<std::string, std::size_t>& variable_map);

    /// @brief Create a process set calculator for a given set of processes with external models
    /// @param processes Processes to create calculator for
    /// @param variable_map A mapping of species names to concentration index
    /// @param external_models External models to include
    CudaProcessSet(
        const std::vector<Process>& processes,
        const std::unordered_map<std::string, std::size_t>& variable_map,
        const std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>>& external_models);

    ~CudaProcessSet()
    {
      micm::cuda::FreeConstData(this->devstruct_);
    }

    /// @brief Set the indexes for the elements of Jacobian matrix before we could copy it to the device;
    /// @brief this will override the "SetJacobianFlatIds" function from the "ProcessSet" class
    /// @param matrix
    void SetJacobianFlatIds(const SparseMatrixPolicy& matrix);

    void AddForcingTerms(const auto& state, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing)
    const requires(VectorizableDense<DenseMatrixPolicy>);

    void SubtractJacobianTerms(
      const auto& state,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
      requires(VectorizableDense<DenseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>);
  };

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::unordered_map<std::string, std::size_t>& variable_map)
      : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(processes, variable_map)
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

    hoststruct.number_of_reactants_size_ = this->number_of_reactants_.size();
    hoststruct.reactant_ids_size_ = this->reactant_ids_.size();
    hoststruct.number_of_products_size_ = this->number_of_products_.size();
    hoststruct.product_ids_size_ = this->product_ids_.size();
    hoststruct.yields_size_ = this->yields_.size();

    // Copy the data from host struct to device struct
    this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::unordered_map<std::string, std::size_t>& variable_map,
      const std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>>& external_models)
      : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(processes, variable_map, external_models)
  {
    // External models are not yet supported on GPU
    if (!external_models.empty())
    {
      throw std::runtime_error("CudaProcessSet does not currently support external models.");
    }

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

    hoststruct.number_of_reactants_size_ = this->number_of_reactants_.size();
    hoststruct.reactant_ids_size_ = this->reactant_ids_.size();
    hoststruct.number_of_products_size_ = this->number_of_products_.size();
    hoststruct.product_ids_size_ = this->product_ids_.size();
    hoststruct.yields_size_ = this->yields_.size();

    // Copy the data from host struct to device struct
    this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline void CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetJacobianFlatIds(const SparseMatrixPolicy& matrix)
  {
    /// This function sets the "jacobian_flat_ids_" member after the structure of Jacobian matrix is known
    ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetJacobianFlatIds(matrix);

    ProcessSetParam hoststruct;
    std::vector<ProcessInfoParam> jacobian_process_info(this->jacobian_process_info_.size());
    std::size_t i_process = 0;
    for (const auto& process_info : this->jacobian_process_info_)
    {
      jacobian_process_info[i_process].process_id_ = process_info.process_id_;
      jacobian_process_info[i_process].independent_id_ = process_info.independent_id_;
      jacobian_process_info[i_process].number_of_dependent_reactants_ = process_info.number_of_dependent_reactants_;
      jacobian_process_info[i_process].number_of_products_ = process_info.number_of_products_;
      ++i_process;
    }
    hoststruct.jacobian_process_info_ = jacobian_process_info.data();
    hoststruct.jacobian_process_info_size_ = jacobian_process_info.size();
    hoststruct.jacobian_reactant_ids_ = this->jacobian_reactant_ids_.data();
    hoststruct.jacobian_reactant_ids_size_ = this->jacobian_reactant_ids_.size();
    hoststruct.jacobian_product_ids_ = this->jacobian_product_ids_.data();
    hoststruct.jacobian_product_ids_size_ = this->jacobian_product_ids_.size();
    hoststruct.jacobian_yields_ = this->jacobian_yields_.data();
    hoststruct.jacobian_yields_size_ = this->jacobian_yields_.size();
    hoststruct.jacobian_flat_ids_ = this->jacobian_flat_ids_.data();
    hoststruct.jacobian_flat_ids_size_ = this->jacobian_flat_ids_.size();

    // Copy the data from host struct to device struct
    micm::cuda::CopyJacobianParams(hoststruct, this->devstruct_);
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline void CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::AddForcingTerms(
    const auto& state,
    const DenseMatrixPolicy& state_variables,
    DenseMatrixPolicy& forcing) const
    requires(VectorizableDense<DenseMatrixPolicy>)
  {
    auto forcing_param = forcing.AsDeviceParam();  // we need to update forcing so it can't be constant and must be an lvalue
    micm::cuda::AddForcingTermsKernelDriver(
        state.rate_constants_.AsDeviceParam(), state_variables.AsDeviceParam(), forcing_param, this->devstruct_);
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline void CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SubtractJacobianTerms(
      const auto& state,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
    requires(VectorizableDense<DenseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  {
    auto jacobian_param =
        jacobian.AsDeviceParam();  // we need to update jacobian so it can't be constant and must be an lvalue
    micm::cuda::SubtractJacobianTermsKernelDriver(
        state.rate_constants_.AsDeviceParam(), state_variables.AsDeviceParam(), jacobian_param, this->devstruct_);
  }
}  // namespace micm
