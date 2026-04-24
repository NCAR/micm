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

    // Parameterized-multiplier rc_index array (static, built once)
    std::size_t* d_mult_rc_indices_ = nullptr;

    // Per-step multiplier values buffer (interleaved layout, grows if needed)
    double*     d_mult_vals_          = nullptr;
    std::size_t d_mult_vals_capacity_ = 0;

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
      FreeDevice(d_mult_rc_indices_);
      FreeDevice(d_mult_vals_);
      d_mult_vals_capacity_ = 0;
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
          d_mult_rc_indices_(std::exchange(other.d_mult_rc_indices_, nullptr)),
          d_mult_vals_(std::exchange(other.d_mult_vals_, nullptr)),
          d_mult_vals_capacity_(std::exchange(other.d_mult_vals_capacity_, 0)),
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

      // Upload parameterized-multiplier rc_indices (static per solver build)
      const auto& mults = cpu_store.parameterized_multipliers_;
      FreeDevice(d_mult_rc_indices_);
      param_.n_multipliers_     = mults.size();
      param_.d_mult_rc_indices_ = nullptr;
      if (!mults.empty())
      {
        std::vector<std::size_t> rc_indices(mults.size());
        for (std::size_t i = 0; i < mults.size(); ++i)
          rc_indices[i] = mults[i].rc_index;
        auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
        CHECK_CUDA_ERROR(
            cudaMallocAsync(&d_mult_rc_indices_, sizeof(std::size_t) * mults.size(), stream), "cudaMalloc");
        CHECK_CUDA_ERROR(
            cudaMemcpyAsync(
                d_mult_rc_indices_, rc_indices.data(), sizeof(std::size_t) * mults.size(), cudaMemcpyHostToDevice, stream),
            "cudaMemcpy");
        param_.d_mult_rc_indices_ = d_mult_rc_indices_;
      }
    }

    /// @brief Evaluate parameterized multipliers on CPU, pack into interleaved layout, and upload.
    ///        Layout: [group * n_mults * L + mult * L + lane]
    /// @return Device pointer to multiplier values, or nullptr if there are no multipliers.
    const double* UploadMultiplierValues(
        const ReactionRateConstantStore& cpu_store,
        const std::vector<Conditions>&   conditions,
        std::size_t                      L)
    {
      const auto& mults = cpu_store.parameterized_multipliers_;
      if (mults.empty())
        return nullptr;

      const std::size_t n_mults  = mults.size();
      const std::size_t n_cells  = conditions.size();
      const std::size_t n_groups = (n_cells + L - 1) / L;
      const std::size_t n_vals   = n_groups * n_mults * L;

      std::vector<double> host_vals(n_vals, 0.0);
      for (std::size_t g = 0; g < n_groups; ++g)
        for (std::size_t i = 0; i < n_mults; ++i)
          for (std::size_t j = 0; j < L; ++j)
          {
            const std::size_t cell = g * L + j;
            if (cell < n_cells)
              host_vals[g * n_mults * L + i * L + j] = mults[i].evaluate(conditions[cell]);
          }

      auto stream = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);
      if (n_vals > d_mult_vals_capacity_)
      {
        FreeDevice(d_mult_vals_);
        CHECK_CUDA_ERROR(cudaMallocAsync(&d_mult_vals_, sizeof(double) * n_vals, stream), "cudaMalloc");
        d_mult_vals_capacity_ = n_vals;
      }
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(d_mult_vals_, host_vals.data(), sizeof(double) * n_vals, cudaMemcpyHostToDevice, stream),
          "cudaMemcpy");
      return d_mult_vals_;
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
