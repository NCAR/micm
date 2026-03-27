// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/process/cuda_process.cuh>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/process/chemical_reaction.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/cuda/solver/cuda_state.hpp>
#include <micm/solver/state.hpp>

#include <cassert>
#include <stdexcept>
#include <vector>

namespace micm
{
  /// @brief A GPU-based implementation for calculating rate constants
  /// @tparam DenseMatrixPolicy Policy for dense matrices (must satisfy CudaMatrix and VectorizableDense concepts)
  template<typename DenseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
  class CudaProcess
  {
   public:
    ProcessParam devstruct_;
    std::vector<CudaRateConstantData> h_rate_constants_;

    CudaProcess() = default;

    CudaProcess(const CudaProcess&) = delete;
    CudaProcess& operator=(const CudaProcess&) = delete;

    CudaProcess(CudaProcess&& other)
        : h_rate_constants_(std::move(other.h_rate_constants_))
    {
      std::swap(this->devstruct_, other.devstruct_);
    }

    CudaProcess& operator=(CudaProcess&& other)
    {
      if (this != &other)
      {
        if (devstruct_.rate_constants_ != nullptr)
          micm::cuda::FreeProcessConstData(devstruct_);
        h_rate_constants_ = std::move(other.h_rate_constants_);
        devstruct_ = {};
        std::swap(this->devstruct_, other.devstruct_);
      }
      return *this;
    }

    ~CudaProcess()
    {
      micm::cuda::FreeProcessConstData(devstruct_);
    }

    /// @brief Construct from a vector of Process objects
    /// @param processes The processes to build rate constant descriptors for
    CudaProcess(const std::vector<Process>& processes);

    /// @brief Calculate rate constants on GPU
    /// @param processes The processes (needed for parameterized reactant evaluation on CPU)
    /// @param state The CUDA solver state (conditions must be synced to device beforehand via SyncConditionsToDevice())
    template<class SparseMatrixPolicy, class LuDecompositionPolicy>
    void CalculateRateConstants(
        const std::vector<Process>& processes,
        CudaState<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>& state) const;

   private:
    /// @brief Pack a rate constant into a CudaRateConstantData descriptor
    /// @param rate_constant The polymorphic rate constant to pack
    /// @return The GPU-portable rate constant descriptor
    static CudaRateConstantData PackRateConstant(const RateConstant* rate_constant);
  };

  // ==========================================================================
  // Implementation
  // ==========================================================================

  template<typename DenseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
  inline CudaProcess<DenseMatrixPolicy>::CudaProcess(const std::vector<Process>& processes)
  {
    // Extract ChemicalReaction objects and pack their rate constants
    for (const auto& process : processes)
    {
      if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
      {
        h_rate_constants_.push_back(PackRateConstant(reaction->rate_constant_.get()));
      }
    }

    // Copy to device
    ProcessParam hoststruct;
    hoststruct.rate_constants_ = h_rate_constants_.data();
    hoststruct.num_reactions_ = h_rate_constants_.size();

    devstruct_ = micm::cuda::CopyProcessConstData(hoststruct);
  }

  template<typename DenseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
  inline CudaRateConstantData CudaProcess<DenseMatrixPolicy>::PackRateConstant(const RateConstant* rate_constant)
  {
    CudaRateConstantData data{};
    data.num_custom_params_ = rate_constant->SizeCustomParameters();

    if (auto* rc = dynamic_cast<const ArrheniusRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Arrhenius;
      data.params_[0] = rc->parameters_.A_;
      data.params_[1] = rc->parameters_.B_;
      data.params_[2] = rc->parameters_.C_;
      data.params_[3] = rc->parameters_.D_;
      data.params_[4] = rc->parameters_.E_;
    }
    else if (auto* rc = dynamic_cast<const TroeRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Troe;
      data.params_[0] = rc->parameters_.k0_A_;
      data.params_[1] = rc->parameters_.k0_B_;
      data.params_[2] = rc->parameters_.k0_C_;
      data.params_[3] = rc->parameters_.kinf_A_;
      data.params_[4] = rc->parameters_.kinf_B_;
      data.params_[5] = rc->parameters_.kinf_C_;
      data.params_[6] = rc->parameters_.Fc_;
      data.params_[7] = rc->parameters_.N_;
    }
    else if (auto* rc = dynamic_cast<const TunnelingRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Tunneling;
      data.params_[0] = rc->parameters_.A_;
      data.params_[1] = rc->parameters_.B_;
      data.params_[2] = rc->parameters_.C_;
    }
    else if (auto* rc = dynamic_cast<const BranchedRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Branched;
      data.params_[0] = rc->parameters_.X_;
      data.params_[1] = rc->parameters_.Y_;
      data.params_[2] = rc->k0_;  // precomputed
      data.params_[3] = rc->z_;   // precomputed
      data.params_[4] =
          (rc->parameters_.branch_ == BranchedRateConstantParameters::Branch::Alkoxy) ? 0.0 : 1.0;
    }
    else if (auto* rc = dynamic_cast<const TernaryChemicalActivationRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::TernaryChemicalActivation;
      data.params_[0] = rc->parameters_.k0_A_;
      data.params_[1] = rc->parameters_.k0_B_;
      data.params_[2] = rc->parameters_.k0_C_;
      data.params_[3] = rc->parameters_.kinf_A_;
      data.params_[4] = rc->parameters_.kinf_B_;
      data.params_[5] = rc->parameters_.kinf_C_;
      data.params_[6] = rc->parameters_.Fc_;
      data.params_[7] = rc->parameters_.N_;
    }
    else if (auto* rc = dynamic_cast<const ReversibleRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Reversible;
      data.params_[0] = rc->parameters_.A_;
      data.params_[1] = rc->parameters_.C_;
      data.params_[2] = rc->parameters_.k_r_;
    }
    else if (auto* rc = dynamic_cast<const SurfaceRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::Surface;
      data.params_[0] = rc->diffusion_coefficient_;
      data.params_[1] = rc->mean_free_speed_factor_;
      data.params_[2] = rc->parameters_.reaction_probability_;
    }
    else if (auto* rc = dynamic_cast<const UserDefinedRateConstant*>(rate_constant))
    {
      data.type_ = CudaRateConstantType::UserDefined;
      data.params_[0] = rc->parameters_.scaling_factor_;
    }
    else
    {
      throw std::runtime_error(
          "CudaProcess: Unsupported rate constant type. "
          "Lambda and TaylorSeries rate constants cannot run on GPU.");
    }

    return data;
  }

  template<typename DenseMatrixPolicy>
    requires(CudaMatrix<DenseMatrixPolicy> && VectorizableDense<DenseMatrixPolicy>)
  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void CudaProcess<DenseMatrixPolicy>::CalculateRateConstants(
      const std::vector<Process>& processes,
      CudaState<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy>& state) const
  {
    const std::size_t num_cells = state.conditions_.size();
    const std::size_t num_reactions = h_rate_constants_.size();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();

    // Verify that the processes vector matches what was used at construction
    std::size_t num_chemical_reactions = 0;
    for (const auto& process : processes)
      if (std::holds_alternative<ChemicalReaction>(process.process_))
        ++num_chemical_reactions;
    assert(num_chemical_reactions == num_reactions && "processes must match the vector used to construct CudaProcess");

    if (num_reactions == 0)
      return;

    auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

    // 1. Conditions (temperature, pressure, air_density) are already on device
    //    via state.SyncConditionsToDevice() called by the caller

    // 2. Pre-compute parameterized reactant factors on CPU
    //    Layout matches rate_constants_ vectorized layout: groups of L cells x num_reactions columns
    std::vector<double> h_fixed_reactants(state.conditions_param_.fixed_reactants_size_, 1.0);

    // Collect chemical reactions
    std::vector<const ChemicalReaction*> reactions;
    for (const auto& process : processes)
    {
      if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
      {
        reactions.push_back(reaction);
      }
    }

    // Compute fixed_reactants in vectorized layout
    for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
    {
      std::size_t offset_rc = i_group * state.rate_constants_.GroupSize();
      std::size_t rate_const_size = std::min(L, state.rate_constants_.NumRows() - (i_group * L));
      for (std::size_t i_rxn = 0; i_rxn < num_reactions; ++i_rxn)
      {
        for (std::size_t i_cell = 0; i_cell < rate_const_size; ++i_cell)
        {
          double fixed = 1.0;
          for (auto& reactant : reactions[i_rxn]->reactants_)
          {
            if (reactant.IsParameterized())
              fixed *= reactant.parameterize_(state.conditions_[i_group * L + i_cell]);
          }
          h_fixed_reactants[offset_rc + i_cell] = fixed;
        }
        offset_rc += L;
      }
    }

    // Copy fixed_reactants to pre-allocated device buffer
    CHECK_CUDA_ERROR(
        cudaMemcpyAsync(
            state.conditions_param_.d_fixed_reactants_,
            h_fixed_reactants.data(),
            sizeof(double) * h_fixed_reactants.size(),
            cudaMemcpyHostToDevice,
            cuda_stream_id),
        "cudaMemcpy");

    // 3. Ensure custom_rate_parameters are on device
    state.custom_rate_parameters_.CopyToDevice();

    // 4. Launch kernel
    auto rate_constants_param = state.rate_constants_.AsDeviceParam();
    auto custom_rate_params = state.custom_rate_parameters_.AsDeviceParam();

    micm::cuda::CalculateRateConstantsKernelDriver(
        state.conditions_param_.d_temperature_,
        state.conditions_param_.d_pressure_,
        state.conditions_param_.d_air_density_,
        custom_rate_params,
        state.conditions_param_.d_fixed_reactants_,
        rate_constants_param,
        devstruct_);
  }

}  // namespace micm
