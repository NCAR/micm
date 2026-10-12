// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

namespace micm
{
  struct TernaryChemicalActivationRateConstantParameters
  {
    /// @brief low-pressure pre-exponential factor
    Real k0_A_ = 1.0;
    /// @brief low-pressure temperature-scaling parameter
    Real k0_B_ = 0.0;
    /// @brief low-pressure exponential factor
    Real k0_C_ = 0.0;
    /// @brief high-pressure pre-exponential factor
    Real kinf_A_ = 1.0;
    /// @brief high-pressure temperature-scaling parameter
    Real kinf_B_ = 0.0;
    /// @brief high-pressure exponential factor
    Real kinf_C_ = 0.0;
    /// @brief TernaryChemicalActivation F_c parameter
    Real Fc_ = 0.6;
    /// @brief TernaryChemicalActivation N parameter
    Real N_ = 1.0;
  };

  /// @brief JPL-19 ternary chemical activation parameters.
  ///        Each of k_0, k_inf and k_int is an Arrhenius-like rate constant:
  ///          k = A * exp(C/T) * (T/D)^B
  ///        Mapping from JPL-19 (Section 2.1):
  ///          k_0:   A = k_0^298,   B = -n, C = 0,  D = 298
  ///          k_inf: A = k_inf^298, B = -m, C = 0,  D = 298
  ///          k_int: A = A,         B = 0,  C = -B, D = unused
  ///        The total rate constant is
  ///          k_total = k_f + k_CA
  ///          k_f     = Troe(k_0, k_inf, [M], Fc, N)
  ///          k_CA    = k_int * (1 - k_f / k_inf)
  struct TernaryChemicalActivationJPL19Parameters
  {
    /// @brief low-pressure limit pre-exponential factor
    Real k0_A_ = 1.0;
    /// @brief low-pressure limit temperature-scaling parameter
    Real k0_B_ = 0.0;
    /// @brief low-pressure limit exponential factor
    Real k0_C_ = 0.0;
    /// @brief low-pressure limit reference temperature
    Real k0_D_ = 298.0;
    /// @brief high-pressure limit pre-exponential factor
    Real kinf_A_ = 1.0;
    /// @brief high-pressure limit temperature-scaling parameter
    Real kinf_B_ = 0.0;
    /// @brief high-pressure limit exponential factor
    Real kinf_C_ = 0.0;
    /// @brief high-pressure limit reference temperature
    Real kinf_D_ = 298.0;
    /// @brief chemical activation (k_int) pre-exponential factor
    Real kint_A_ = 1.0;
    /// @brief chemical activation (k_int) temperature-scaling parameter
    Real kint_B_ = 0.0;
    /// @brief chemical activation (k_int) exponential factor
    Real kint_C_ = 0.0;
    /// @brief chemical activation (k_int) reference temperature (does not affect the result when kint_B_ = 0)
    Real kint_D_ = 298.0;
    /// @brief F_c parameter
    Real Fc_ = 0.6;
    /// @brief N parameter
    Real N_ = 1.0;
  };
}  // namespace micm
