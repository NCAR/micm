// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file cuda_rate_constant_kernel.cu
/// @brief GPU kernel for analytic rate constant calculation.
///
/// One thread per grid cell.  Writes directly into the VectorMatrix interleaved
/// layout of rate_constants_.  Lambda entries (written by EvaluateCpuRateConstants on the
/// CPU) are not touched.

#define _USE_MATH_DEFINES

#include <micm/cuda/process/cuda_rate_constant_kernel.cuh>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/util/types.hpp>

#include <cmath>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

namespace micm::cuda
{
  /// @brief Compute all analytic rate constants for one thread and write to rc_base (interleaved).
  /// @param cp_base   Group base pointer into custom_rate_parameters_ for this thread
  /// @param rc_base   Group base pointer into rate_constants_ for this thread
  /// @param mv_base   Group base pointer into mult_vals buffer [mult * L + lane]; nullptr if no multipliers
  /// @param local_tid Lane index within the VectorMatrix group (tid % L)
  __device__ __forceinline__ static void CalculateRatesForThread(
      const CudaReactionRateStoreParam& store,
      Real temperature,
      Real pressure,
      Real air_density,
      const Real* cp_base,
      Real* rc_base,
      const Real* mv_base,
      Index local_tid,
      Index L)
  {
    auto out = [&](Index offset, Index i) -> Real* { return &rc_base[(offset + i) * L + local_tid]; };

    for (Index i = 0; i < store.n_arrhenius_; ++i)
    {
      *out(0, i) = micm::CalculateArrhenius(store.d_arrhenius_[i], temperature, pressure);
    }
    for (Index i = 0; i < store.n_troe_; ++i)
    {
      *out(store.troe_offset_, i) = micm::CalculateTroe(store.d_troe_[i], temperature, air_density);
    }
    for (Index i = 0; i < store.n_ternary_; ++i)
    {
      *out(store.ternary_offset_, i) =
          micm::CalculateTernaryChemicalActivation(store.d_ternary_[i], temperature, air_density);
    }
    for (Index i = 0; i < store.n_branched_; ++i)
    {
      *out(store.branched_offset_, i) = micm::CalculateBranched(store.d_branched_[i], temperature, air_density);
    }
    for (Index i = 0; i < store.n_tunneling_; ++i)
    {
      *out(store.tunneling_offset_, i) = micm::CalculateTunneling(store.d_tunneling_[i], temperature);
    }
    for (Index i = 0; i < store.n_taylor_; ++i)
    {
      *out(store.taylor_offset_, i) = micm::CalculateTaylorSeries(store.d_taylor_[i], temperature, pressure);
    }
    for (Index i = 0; i < store.n_reversible_; ++i)
    {
      *out(store.reversible_offset_, i) = micm::CalculateReversible(store.d_reversible_[i], temperature);
    }
    for (Index i = 0; i < store.n_user_defined_; ++i)
    {
      const micm::UserDefinedRateConstantData& p = store.d_user_defined_[i];
      *out(store.user_defined_offset_, i) = micm::CalculateUserDefined(p, cp_base[p.custom_param_index_ * L + local_tid]);
    }
    for (Index i = 0; i < store.n_surface_; ++i)
    {
      const micm::SurfaceRateConstantData& p = store.d_surface_[i];
      *out(store.surface_offset_, i) = micm::CalculateSurfaceOne(
          p,
          temperature,
          cp_base[p.custom_param_base_index_ * L + local_tid],
          cp_base[(p.custom_param_base_index_ + 1) * L + local_tid]);
    }
    for (Index i = 0; i < store.n_multipliers_; ++i)
    {
      rc_base[store.d_mult_rc_indices_[i] * L + local_tid] *= mv_base[i * L + local_tid];
    }
  }

  /// @brief Calculate all analytic rate constants for every grid cell (one thread per cell).
  ///
  /// Interleaved layout: element (cell tid, col k) = base[k * L + local_tid]
  ///   where base = group_id * group_size, local_tid = tid % L.
  __global__ void CalculateRateConstantsKernel(
      const CudaReactionRateStoreParam store,
      const micm::Conditions* d_conditions,
      Real* d_rc,
      Index rc_group_size,
      const Real* d_cp,
      Index cp_group_size,
      const Real* d_mult_vals,
      Index mult_group_size,
      Index n_cells,
      Index L)
  {
    const Index tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    const Index group_id = tid / L;
    const Index local_tid = tid % L;

    if (tid >= n_cells)
    {
      return;
    }

    const Real temperature = d_conditions[tid].temperature_;
    const Real pressure = d_conditions[tid].pressure_;
    const Real air_density = d_conditions[tid].air_density_;

    CalculateRatesForThread(
        store,
        temperature,
        pressure,
        air_density,
        d_cp + group_id * cp_group_size,
        d_rc + group_id * rc_group_size,
        d_mult_vals + group_id * mult_group_size,
        local_tid,
        L);
  }

  void CalculateRateConstantsKernelDriver(
      const CudaReactionRateStoreParam& store_param,
      const micm::Conditions* d_conditions,
      CudaMatrixParam& rc_param,
      const CudaMatrixParam& cp_param,
      const Real* d_mult_vals)
  {
    const Index n_cells = rc_param.number_of_grid_cells_;
    if (n_cells == 0)
    {
      return;
    }

    const Index L = rc_param.vector_length_;
    const Index n_groups = (n_cells + L - 1) / L;
    const Index rc_group_size = rc_param.number_of_elements_ / n_groups;
    const Index cp_group_size = (cp_param.number_of_elements_ > 0) ? cp_param.number_of_elements_ / n_groups : 0;
    const Index mult_group_size = store_param.n_multipliers_ * L;

    const Index number_of_blocks = (n_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;

    CalculateRateConstantsKernel<<<
        number_of_blocks,
        BLOCK_SIZE,
        0,
        micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
        store_param,
        d_conditions,
        rc_param.d_data_,
        rc_group_size,
        cp_param.d_data_,
        cp_group_size,
        d_mult_vals,
        mult_group_size,
        n_cells,
        L);
  }
}  // namespace micm::cuda
