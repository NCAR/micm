// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file cuda_rate_constant_kernel.cu
/// @brief GPU kernel for two-phase analytic rate constant calculation (Phase 2).
///
/// Each CUDA thread handles one grid cell.  For every analytic reaction type the
/// thread calls the corresponding MICM_HOST_DEVICE batch function with n=1 and
/// writes the result directly into the VectorMatrix interleaved layout of
/// rate_constants_.
///
/// Lambda entries (written on the CPU by EvaluateCpuRates) occupy positions
/// beyond surface_offset in rate_constants_ and are NOT touched by this kernel.
///
/// UserDefined and Surface rate constants read their inputs from the
/// interleaved custom_rate_parameters_ layout, so they must be accessed via
/// the stride formulae rather than through a gathered temporary buffer.

#define _USE_MATH_DEFINES

#include <micm/cuda/process/cuda_rate_constant_kernel.cuh>
#include <micm/cuda/util/cuda_util.cuh>
#include <micm/process/rate_constant/rate_constant_functions.hpp>

#include <math.h>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace micm
{
  namespace cuda
  {
    /// @brief GPU kernel: calculate all analytic rate constants for every grid cell.
    ///
    /// Thread mapping:
    ///   tid       = global thread index = cell index
    ///   group_id  = tid / L (VectorMatrix group)
    ///   local_tid = tid % L (lane within group)
    ///
    /// Layout in rate_constants_ (n_cells rows × n_rc columns, interleaved):
    ///   element (cell tid, reaction k) = d_rc[rc_group_base + k * L + local_tid]
    ///   where rc_group_base = group_id * (total_elements / n_groups)
    ///
    /// Layout in custom_rate_parameters_ (n_cells rows × n_cp columns, interleaved):
    ///   element (cell tid, param j)    = d_cp[cp_group_base + j * L + local_tid]
    __global__ void CalculateRateConstantsKernel(
        const CudaReactionRateStoreParam store,
        const micm::Conditions*          d_conditions,
        double*                          d_rc,
        std::size_t                      rc_group_size,
        const double*                    d_cp,
        std::size_t                      cp_group_size,
        std::size_t                      n_cells,
        std::size_t                      L)
    {
      const std::size_t tid       = blockIdx.x * BLOCK_SIZE + threadIdx.x;
      const std::size_t group_id  = tid / L;
      const std::size_t local_tid = tid % L;

      if (tid >= n_cells)
        return;

      const double temperature = d_conditions[tid].temperature_;
      const double pressure    = d_conditions[tid].pressure_;
      const double air_density = d_conditions[tid].air_density_;

      // Pointers to the start of this thread's group in the interleaved arrays
      const double* cp_base = d_cp + group_id * cp_group_size;
      double*       rc_base = d_rc + group_id * rc_group_size;

      // Macro: call a batch function with n=1 for reaction i, write at offset+i
#define WRITE_RC(offset, i, val) rc_base[((offset) + (i)) * L + local_tid] = (val)

      // ----------------------------------------------------------------
      // Arrhenius (offset = 0)
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_arrhenius_; ++i)
      {
        double val;
        micm::CalculateArrhenius(store.d_arrhenius_ + i, 1, temperature, pressure, &val);
        WRITE_RC(0, i, val);
      }

      // ----------------------------------------------------------------
      // Troe
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_troe_; ++i)
      {
        double val;
        micm::CalculateTroe(store.d_troe_ + i, 1, temperature, air_density, &val);
        WRITE_RC(store.troe_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // TernaryChemicalActivation
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_ternary_; ++i)
      {
        double val;
        micm::CalculateTernaryChemicalActivation(store.d_ternary_ + i, 1, temperature, air_density, &val);
        WRITE_RC(store.ternary_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // Branched
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_branched_; ++i)
      {
        double val;
        micm::CalculateBranched(store.d_branched_ + i, 1, temperature, air_density, &val);
        WRITE_RC(store.branched_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // Tunneling
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_tunneling_; ++i)
      {
        double val;
        micm::CalculateTunneling(store.d_tunneling_ + i, 1, temperature, &val);
        WRITE_RC(store.tunneling_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // TaylorSeries
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_taylor_; ++i)
      {
        double val;
        micm::CalculateTaylorSeries(store.d_taylor_ + i, 1, temperature, pressure, &val);
        WRITE_RC(store.taylor_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // Reversible
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_reversible_; ++i)
      {
        double val;
        micm::CalculateReversible(store.d_reversible_ + i, 1, temperature, &val);
        WRITE_RC(store.reversible_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // UserDefined — read from interleaved custom_rate_parameters_
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_user_defined_; ++i)
      {
        const micm::UserDefinedRateConstantData& p   = store.d_user_defined_[i];
        double                                   val = cp_base[p.custom_param_index_ * L + local_tid] * p.scaling_factor_;
        WRITE_RC(store.user_defined_offset_, i, val);
      }

      // ----------------------------------------------------------------
      // Surface — read radius and num_conc from interleaved custom_rate_parameters_
      // ----------------------------------------------------------------
      for (std::size_t i = 0; i < store.n_surface_; ++i)
      {
        const micm::SurfaceRateConstantData& p        = store.d_surface_[i];
        double                               radius   = cp_base[p.custom_param_base_index_ * L + local_tid];
        double                               num_conc = cp_base[(p.custom_param_base_index_ + 1) * L + local_tid];
        double mean_free_speed = sqrt(p.mean_free_speed_factor_ * temperature);
        double val             = 4.0 * num_conc * M_PI * radius * radius /
                     (radius / p.diffusion_coefficient_ + 4.0 / (mean_free_speed * p.reaction_probability_));
        WRITE_RC(store.surface_offset_, i, val);
      }

#undef WRITE_RC
    }

    void CalculateRateConstantsKernelDriver(
        const CudaReactionRateStoreParam& store_param,
        const micm::Conditions*           d_conditions,
        CudaMatrixParam&                  rc_param,
        const CudaMatrixParam&            cp_param)
    {
      const std::size_t n_cells = rc_param.number_of_grid_cells_;
      if (n_cells == 0)
        return;

      const std::size_t L             = rc_param.vector_length_;
      const std::size_t n_groups      = (n_cells + L - 1) / L;
      const std::size_t rc_group_size = rc_param.number_of_elements_ / n_groups;
      const std::size_t cp_group_size = (cp_param.number_of_elements_ > 0) ? cp_param.number_of_elements_ / n_groups : 0;

      const std::size_t number_of_blocks = (n_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;

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
          n_cells,
          L);
    }
  }  // namespace cuda
}  // namespace micm
