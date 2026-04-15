// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/process/cuda_process_set.cuh>
#include <micm/cuda/process/cuda_rate_constant_kernel.cuh>
#include <micm/cuda/process/cuda_reaction_rate_store.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/process/process_set.hpp>
#include <micm/process/reaction_rate_store.hpp>

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

    /// GPU-resident analytic rate constant parameter store (built once per solver build)
    CudaReactionRateStore cuda_rate_store_;

    CudaProcessSet() = default;

    CudaProcessSet(const CudaProcessSet&) = delete;
    CudaProcessSet& operator=(const CudaProcessSet&) = delete;
    CudaProcessSet(CudaProcessSet&& other)
        : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(std::move(other)),
          cuda_rate_store_(std::move(other.cuda_rate_store_))
    {
      std::swap(this->devstruct_, other.devstruct_);
    };
    CudaProcessSet& operator=(CudaProcessSet&& other)
    {
      ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::operator=(std::move(other));
      std::swap(this->devstruct_, other.devstruct_);
      cuda_rate_store_ = std::move(other.cuda_rate_store_);
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

    /// @brief Upload all analytic parameter arrays from cpu_store to device memory.
    ///        Called once by Solver after ReactionRateStore is built.
    void BuildCudaStore(const ReactionRateStore& cpu_store)
    {
      cuda_rate_store_.BuildFrom(cpu_store);
    }

    /// @brief GPU-accelerated rate constant calculation.
    ///
    ///   1. Evaluate any lambda entries on CPU; upload rate_constants_ to device.
    ///   2. Upload conditions and custom_rate_parameters_ to device.
    ///   3. Launch CalculateRateConstantsKernel to fill analytic slots on device.
    ///   4. If parameterized multipliers exist: download, apply on CPU, re-upload.
    ///
    /// After this call, device rate_constants_ is fully populated for the current step.
    template<class StatePolicy>
    void GpuCalculateRateConstants(const ReactionRateStore& cpu_store, StatePolicy& state)
      requires(CudaMatrix<typename StatePolicy::DenseMatrixPolicyType> &&
               VectorizableDense<typename StatePolicy::DenseMatrixPolicyType>)
    {
      // CPU lambda evaluation
      if (!cpu_store.lambda_entries.empty())
      {
        ReactionRateStore::EvaluateCpuRates(cpu_store, state);
        // Upload lambda values (analytic slots carry stale data; kernel overwrites them)
        state.rate_constants_.CopyToDevice();
      }

      // Upload per-step transient data
      const Conditions* d_conditions = cuda_rate_store_.UploadConditions(state.conditions_);
      state.custom_rate_parameters_.CopyToDevice();

      // GPU analytic calculation
      auto rc_param = state.rate_constants_.AsDeviceParam();
      auto cp_param = state.custom_rate_parameters_.AsDeviceParam();
      micm::cuda::CalculateRateConstantsKernelDriver(cuda_rate_store_.GetParam(), d_conditions, rc_param, cp_param);

      // Parameterized multipliers are CPU-only (std::function); apply via round-trip if needed
      if (!cpu_store.parameterized_multipliers.empty())
      {
        state.rate_constants_.CopyToHost();
        const std::size_t n_cells = state.rate_constants_.NumRows();
        for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
        {
          const auto& cond = state.conditions_[i_cell];
          for (const auto& mult : cpu_store.parameterized_multipliers)
            state.rate_constants_[i_cell][mult.rc_index] *= mult.evaluate(cond);
        }
        state.rate_constants_.CopyToDevice();
      }
    }

    /// @brief Set the indexes for the elements of Jacobian matrix before we could copy it to the device;
    /// @brief this will override the "SetJacobianFlatIds" function from the "ProcessSet" class
    /// @param matrix
    void SetJacobianFlatIds(const SparseMatrixPolicy& matrix);

    /// @brief Marks species rows that should be treated as algebraic (constraints replace ODE rows).
    ///        Updates algebraic variable IDs after `ProcessSetParam` construction.
    ///        If algebraic variable IDs are not set post-construction, then this function may not be
    ///        necessary.
    /// @param variable_ids Set of variable ids whose forcing/Jacobian rows should not receive kinetic contributions
    void SetAlgebraicVariableIds(const std::set<std::size_t>& variable_ids);

    void AddForcingTerms(const auto& state, const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& forcing) const
      requires(VectorizableDense<DenseMatrixPolicy>);

    void SubtractJacobianTerms(const auto& state, const DenseMatrixPolicy& state_variables, SparseMatrixPolicy& jacobian)
        const
      requires(VectorizableDense<DenseMatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>);

   private:
    void InitDevStruct();
  };

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline void CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::InitDevStruct()
  {
    ProcessSetParam hoststruct;
    hoststruct.number_of_reactants_      = this->number_of_reactants_.data();
    hoststruct.reactant_ids_             = this->reactant_ids_.data();
    hoststruct.number_of_products_       = this->number_of_products_.data();
    hoststruct.product_ids_              = this->product_ids_.data();
    hoststruct.yields_                   = this->yields_.data();
    hoststruct.is_algebraic_variable_    = this->is_algebraic_variable_.data();
    hoststruct.number_of_reactants_size_ = this->number_of_reactants_.size();
    hoststruct.reactant_ids_size_        = this->reactant_ids_.size();
    hoststruct.number_of_products_size_  = this->number_of_products_.size();
    hoststruct.product_ids_size_         = this->product_ids_.size();
    hoststruct.yields_size_              = this->yields_.size();
    hoststruct.algebraic_variable_size_  = this->is_algebraic_variable_.size();
    this->devstruct_ = micm::cuda::CopyConstData(hoststruct);
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::unordered_map<std::string, std::size_t>& variable_map)
      : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(processes, variable_map)
  {
    InitDevStruct();
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && CudaMatrix<SparseMatrixPolicy>)
  inline CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::CudaProcessSet(
      const std::vector<Process>& processes,
      const std::unordered_map<std::string, std::size_t>& variable_map,
      const std::vector<ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>>& external_models)
      : ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>(processes, variable_map, external_models)
  {
    if (!external_models.empty())
      throw std::runtime_error("CudaProcessSet does not currently support external models.");
    InitDevStruct();
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
  inline void CudaProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetAlgebraicVariableIds(
      const std::set<std::size_t>& variable_ids)
  {
    // Update the host-side is_algebraic_variable_ array
    ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>::SetAlgebraicVariableIds(variable_ids);

    // Then update the device memory
    ProcessSetParam hoststruct;
    hoststruct.is_algebraic_variable_ = this->is_algebraic_variable_.data();
    hoststruct.algebraic_variable_size_ = this->is_algebraic_variable_.size();

    // Copy the data from host struct to device struct
    micm::cuda::CopyAlgebraicVariableParams(hoststruct, this->devstruct_);
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
