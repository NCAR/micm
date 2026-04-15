// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/taylor_series_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/system/conditions.hpp>

#include <cstddef>

/// @brief GPU-side mirror of ReactionRateStore analytic parameters.
///
/// All pointer members must be device pointers when passed to
/// CalculateRateConstantsKernelDriver.  The offset fields are pre-computed
/// cumulative sums matching the ReactionRateStore offset helpers.
struct CudaReactionRateStoreParam
{
  // Analytic parameter arrays (device pointers, read-only during kernel)
  const micm::ArrheniusRateConstantParameters*                 d_arrhenius_    = nullptr;
  const micm::TroeRateConstantParameters*                      d_troe_         = nullptr;
  const micm::TernaryChemicalActivationRateConstantParameters* d_ternary_      = nullptr;
  const micm::BranchedRateConstantParameters*                  d_branched_     = nullptr;
  const micm::TunnelingRateConstantParameters*                 d_tunneling_    = nullptr;
  const micm::TaylorSeriesRateConstantParameters*              d_taylor_       = nullptr;
  const micm::ReversibleRateConstantParameters*                d_reversible_   = nullptr;
  const micm::UserDefinedRateConstantData*                     d_user_defined_ = nullptr;
  const micm::SurfaceRateConstantData*                         d_surface_      = nullptr;

  // Counts for each type
  std::size_t n_arrhenius_    = 0;
  std::size_t n_troe_         = 0;
  std::size_t n_ternary_      = 0;
  std::size_t n_branched_     = 0;
  std::size_t n_tunneling_    = 0;
  std::size_t n_taylor_       = 0;
  std::size_t n_reversible_   = 0;
  std::size_t n_user_defined_ = 0;
  std::size_t n_surface_      = 0;

  // Pre-computed contiguous-block offsets into rate_constants_[cell]
  // (same semantics as the inline helpers in ReactionRateStore)
  std::size_t troe_offset_         = 0;  // = n_arrhenius
  std::size_t ternary_offset_      = 0;
  std::size_t branched_offset_     = 0;
  std::size_t tunneling_offset_    = 0;
  std::size_t taylor_offset_       = 0;
  std::size_t reversible_offset_   = 0;
  std::size_t user_defined_offset_ = 0;
  std::size_t surface_offset_      = 0;
};

namespace micm
{
  namespace cuda
  {
    /// @brief Host-callable driver that launches the rate constant kernel on the GPU.
    ///
    /// Phase 2 of the two-phase calculation: computes all analytic rate constants for
    /// every grid cell in parallel and writes them into the correct type-block offsets
    /// of rate_constants_ (VectorMatrix interleaved layout).
    ///
    /// Lambda entries (written by EvaluateCpuRates / CopyToDevice before this call)
    /// are NOT touched by the kernel — they occupy positions beyond surface_offset.
    ///
    /// @param store_param   Device-side analytic parameter arrays
    /// @param d_conditions  Device array of Conditions structs, length ≥ n_cells
    /// @param rc_param      Rate constants matrix device param (interleaved layout)
    /// @param cp_param      Custom rate parameters matrix device param (interleaved layout)
    void CalculateRateConstantsKernelDriver(
        const CudaReactionRateStoreParam& store_param,
        const micm::Conditions*           d_conditions,
        CudaMatrixParam&                  rc_param,
        const CudaMatrixParam&            cp_param);
  }  // namespace cuda
}  // namespace micm
