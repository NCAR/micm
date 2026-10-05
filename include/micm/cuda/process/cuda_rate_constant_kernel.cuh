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
#include <micm/util/types.hpp>

#include <cstddef>

/// @brief Device-side mirror of ReactionRateConstantStore for use in GPU kernels.
///
/// All pointers must be device pointers.  Offsets match ReactionRateConstantStore's
/// offset helpers (cumulative type counts).
struct CudaReactionRateStoreParam
{
  // Device parameter arrays (read-only in kernel)
  const micm::ArrheniusRateConstantParameters* d_arrhenius_ = nullptr;
  const micm::TroeRateConstantParameters* d_troe_ = nullptr;
  const micm::TernaryChemicalActivationRateConstantParameters* d_ternary_ = nullptr;
  const micm::BranchedRateConstantParameters* d_branched_ = nullptr;
  const micm::TunnelingRateConstantParameters* d_tunneling_ = nullptr;
  const micm::TaylorSeriesRateConstantParameters* d_taylor_ = nullptr;
  const micm::ReversibleRateConstantParameters* d_reversible_ = nullptr;
  const micm::UserDefinedRateConstantData* d_user_defined_ = nullptr;
  const micm::SurfaceRateConstantData* d_surface_ = nullptr;

  // Reaction counts per type
  micm::Index n_arrhenius_ = 0;
  micm::Index n_troe_ = 0;
  micm::Index n_ternary_ = 0;
  micm::Index n_branched_ = 0;
  micm::Index n_tunneling_ = 0;
  micm::Index n_taylor_ = 0;
  micm::Index n_reversible_ = 0;
  micm::Index n_user_defined_ = 0;
  micm::Index n_surface_ = 0;

  // Offsets into rate_constants_[cell] (cumulative type counts)
  micm::Index troe_offset_ = 0;
  micm::Index ternary_offset_ = 0;
  micm::Index branched_offset_ = 0;
  micm::Index tunneling_offset_ = 0;
  micm::Index taylor_offset_ = 0;
  micm::Index reversible_offset_ = 0;
  micm::Index user_defined_offset_ = 0;
  micm::Index surface_offset_ = 0;

  // Parameterized-reactant multipliers (static per solver build)
  const micm::Index* d_mult_rc_indices_ = nullptr;
  micm::Index n_multipliers_ = 0;
};

namespace micm::cuda
{
  /// @brief Launch the rate constant kernel.  Lambda entries are not touched.
  /// @param d_mult_vals  Per-step interleaved multiplier values [group * n_mults * L + mult * L + lane].
  ///                     Nullptr when n_multipliers_ == 0.
  void CalculateRateConstantsKernelDriver(
      const CudaReactionRateStoreParam& store_param,
      const micm::Conditions* d_conditions,
      CudaMatrixParam& rc_param,
      const CudaMatrixParam& cp_param,
      const Real* d_mult_vals);
}  // namespace micm::cuda
