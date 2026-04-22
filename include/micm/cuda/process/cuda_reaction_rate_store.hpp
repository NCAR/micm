// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file cuda_reaction_rate_store.hpp
/// @brief GPU-resident store for analytic rate constant parameters.
///
/// CudaReactionRateStore uploads all analytic parameter arrays from a
/// ReactionRateConstantStore to device memory once at construction (via BuildFrom).
/// Per-step transient data (conditions, custom_params) is uploaded each step
/// inside GpuCalculateRateConstants on CudaProcessSet.

#include <micm/cuda/process/cuda_rate_constant_kernel.cuh>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/process/reaction_rate_store.hpp>
#include <micm/system/conditions.hpp>

namespace micm
{
  /// @brief GPU-resident mirror of ReactionRateConstantStore analytic data.
  ///
  /// Constructed once per solver build; never modified during a run.
  /// The device conditions buffer grows on demand (amortised allocation).
  class CudaReactionRateStore
  {
   private:
    // ----------------------------------------------------------------
    // Owned device pointers — freed in destructor
    // ----------------------------------------------------------------
    ArrheniusRateConstantParameters*                 d_arrhenius_    = nullptr;
    TroeRateConstantParameters*                      d_troe_         = nullptr;
    TernaryChemicalActivationRateConstantParameters* d_ternary_      = nullptr;
    BranchedRateConstantParameters*                  d_branched_     = nullptr;
    TunnelingRateConstantParameters*                 d_tunneling_    = nullptr;
    TaylorSeriesRateConstantParameters*              d_taylor_       = nullptr;
    ReversibleRateConstantParameters*                d_reversible_   = nullptr;
    UserDefinedRateConstantData*                     d_user_defined_ = nullptr;
    SurfaceRateConstantData*                         d_surface_      = nullptr;

    // Device buffer for per-step conditions (grows if needed)
    Conditions* d_conditions_          = nullptr;
    std::size_t d_conditions_capacity_ = 0;

    // Cached kernel param struct (populated by BuildFrom)
    CudaReactionRateStoreParam param_{};

    // ----------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------
    template<class T>
    static void ReallocAndUpload(T*& d_ptr, const std::vector<T>& host_vec)
    {
      FreeDevice(d_ptr);
      if (host_vec.empty())
        return;
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      CHECK_CUDA_ERROR(cudaMallocAsync(&d_ptr, sizeof(T) * host_vec.size(), stream), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(d_ptr, host_vec.data(), sizeof(T) * host_vec.size(), cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
    }

    template<class T>
    static void FreeDevice(T*& d_ptr)
    {
      if (d_ptr)
      {
        auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
        CHECK_CUDA_ERROR(cudaFreeAsync(d_ptr, stream), "cudaFree");
        d_ptr = nullptr;
      }
    }

    void FreeAll()
    {
      FreeDevice(d_arrhenius_);
      FreeDevice(d_troe_);
      FreeDevice(d_ternary_);
      FreeDevice(d_branched_);
      FreeDevice(d_tunneling_);
      FreeDevice(d_taylor_);
      FreeDevice(d_reversible_);
      FreeDevice(d_user_defined_);
      FreeDevice(d_surface_);
      FreeDevice(d_conditions_);
      d_conditions_capacity_ = 0;
    }

   public:
    CudaReactionRateStore() = default;

    CudaReactionRateStore(const CudaReactionRateStore&)            = delete;
    CudaReactionRateStore& operator=(const CudaReactionRateStore&) = delete;

    CudaReactionRateStore(CudaReactionRateStore&& other) noexcept
        : d_arrhenius_(std::exchange(other.d_arrhenius_, nullptr)),
          d_troe_(std::exchange(other.d_troe_, nullptr)),
          d_ternary_(std::exchange(other.d_ternary_, nullptr)),
          d_branched_(std::exchange(other.d_branched_, nullptr)),
          d_tunneling_(std::exchange(other.d_tunneling_, nullptr)),
          d_taylor_(std::exchange(other.d_taylor_, nullptr)),
          d_reversible_(std::exchange(other.d_reversible_, nullptr)),
          d_user_defined_(std::exchange(other.d_user_defined_, nullptr)),
          d_surface_(std::exchange(other.d_surface_, nullptr)),
          d_conditions_(std::exchange(other.d_conditions_, nullptr)),
          d_conditions_capacity_(std::exchange(other.d_conditions_capacity_, 0)),
          param_(std::exchange(other.param_, {}))
    {
    }

    CudaReactionRateStore& operator=(CudaReactionRateStore&& other) noexcept
    {
      if (this != &other)
      {
        FreeAll();
        new (this) CudaReactionRateStore(std::move(other));
      }
      return *this;
    }

    ~CudaReactionRateStore()
    {
      FreeAll();
    }

    /// @brief Upload all analytic parameter arrays from cpu_store to device memory.
    ///
    ///        Called once after the ReactionRateConstantStore is built in Solver's constructor.
    ///        Any previous device allocations are freed before re-uploading.
    void BuildFrom(const ReactionRateConstantStore& cpu_store)
    {
      ReallocAndUpload(d_arrhenius_,    cpu_store.arrhenius_);
      ReallocAndUpload(d_troe_,         cpu_store.troe_);
      ReallocAndUpload(d_ternary_,      cpu_store.ternary_);
      ReallocAndUpload(d_branched_,     cpu_store.branched_);
      ReallocAndUpload(d_tunneling_,    cpu_store.tunneling_);
      ReallocAndUpload(d_taylor_,       cpu_store.taylor_);
      ReallocAndUpload(d_reversible_,   cpu_store.reversible_);
      ReallocAndUpload(d_user_defined_, cpu_store.user_defined_);
      ReallocAndUpload(d_surface_,      cpu_store.surface_);

      // Populate the kernel param struct
      param_.d_arrhenius_    = d_arrhenius_;
      param_.d_troe_         = d_troe_;
      param_.d_ternary_      = d_ternary_;
      param_.d_branched_     = d_branched_;
      param_.d_tunneling_    = d_tunneling_;
      param_.d_taylor_       = d_taylor_;
      param_.d_reversible_   = d_reversible_;
      param_.d_user_defined_ = d_user_defined_;
      param_.d_surface_      = d_surface_;

      param_.n_arrhenius_    = cpu_store.arrhenius_.size();
      param_.n_troe_         = cpu_store.troe_.size();
      param_.n_ternary_      = cpu_store.ternary_.size();
      param_.n_branched_     = cpu_store.branched_.size();
      param_.n_tunneling_    = cpu_store.tunneling_.size();
      param_.n_taylor_       = cpu_store.taylor_.size();
      param_.n_reversible_   = cpu_store.reversible_.size();
      param_.n_user_defined_ = cpu_store.user_defined_.size();
      param_.n_surface_      = cpu_store.surface_.size();

      param_.troe_offset_         = cpu_store.troe_offset();
      param_.ternary_offset_      = cpu_store.ternary_offset();
      param_.branched_offset_     = cpu_store.branched_offset();
      param_.tunneling_offset_    = cpu_store.tunneling_offset();
      param_.taylor_offset_       = cpu_store.taylor_offset();
      param_.reversible_offset_   = cpu_store.reversible_offset();
      param_.user_defined_offset_ = cpu_store.user_defined_offset();
      param_.surface_offset_      = cpu_store.surface_offset();
    }

    /// @brief Upload the current conditions array to device, growing the buffer if needed.
    /// @return Device pointer valid until the next call to UploadConditions.
    const Conditions* UploadConditions(const std::vector<Conditions>& conditions)
    {
      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      if (conditions.size() > d_conditions_capacity_)
      {
        FreeDevice(d_conditions_);
        CHECK_CUDA_ERROR(
            cudaMallocAsync(&d_conditions_, sizeof(Conditions) * conditions.size(), stream), "cudaMalloc");
        d_conditions_capacity_ = conditions.size();
      }
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              d_conditions_,
              conditions.data(),
              sizeof(Conditions) * conditions.size(),
              cudaMemcpyHostToDevice,
              stream),
          "cudaMemcpy");
      return d_conditions_;
    }

    const CudaReactionRateStoreParam& GetParam() const
    {
      return param_;
    }
  };

}  // namespace micm
