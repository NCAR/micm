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
}  // namespace micm
