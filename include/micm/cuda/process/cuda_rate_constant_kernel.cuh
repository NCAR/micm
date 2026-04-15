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

/// @brief Device-side mirror of ReactionRateStore for use in GPU kernels.
///
/// All pointers must be device pointers.  Offsets match ReactionRateStore's
/// offset helpers (cumulative type counts).
struct CudaReactionRateStoreParam
{
  // Device parameter arrays (read-only in kernel)
  const micm::ArrheniusRateConstantParameters*                 d_arrhenius_    = nullptr;
  const micm::TroeRateConstantParameters*                      d_troe_         = nullptr;
  const micm::TernaryChemicalActivationRateConstantParameters* d_ternary_      = nullptr;
  const micm::BranchedRateConstantParameters*                  d_branched_     = nullptr;
  const micm::TunnelingRateConstantParameters*                 d_tunneling_    = nullptr;
  const micm::TaylorSeriesRateConstantParameters*              d_taylor_       = nullptr;
  const micm::ReversibleRateConstantParameters*                d_reversible_   = nullptr;
  const micm::UserDefinedRateConstantData*                     d_user_defined_ = nullptr;
  const micm::SurfaceRateConstantData*                         d_surface_      = nullptr;

  // Reaction counts per type
  std::size_t n_arrhenius_    = 0;
  std::size_t n_troe_         = 0;
  std::size_t n_ternary_      = 0;
  std::size_t n_branched_     = 0;
  std::size_t n_tunneling_    = 0;
  std::size_t n_taylor_       = 0;
  std::size_t n_reversible_   = 0;
  std::size_t n_user_defined_ = 0;
  std::size_t n_surface_      = 0;

  // Offsets into rate_constants_[cell] (cumulative type counts)
  std::size_t troe_offset_         = 0;
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
    /// @brief Launch the rate constant kernel.  Lambda entries are not touched.
    void CalculateRateConstantsKernelDriver(
        const CudaReactionRateStoreParam& store_param,
        const micm::Conditions*           d_conditions,
        CudaMatrixParam&                  rc_param,
        const CudaMatrixParam&            cp_param);
  }  // namespace cuda
}  // namespace micm
