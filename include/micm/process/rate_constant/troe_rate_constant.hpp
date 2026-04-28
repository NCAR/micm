// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace micm
{
  struct TroeRateConstantParameters
  {
    /// @brief low-pressure pre-exponential factor
    double k0_A_ = 1.0;
    /// @brief low-pressure temperature-scaling parameter
    double k0_B_ = 0.0;
    /// @brief low-pressure exponential factor
    double k0_C_ = 0.0;
    /// @brief high-pressure pre-exponential factor
    double kinf_A_ = 1.0;
    /// @brief high-pressure temperature-scaling parameter
    double kinf_B_ = 0.0;
    /// @brief high-pressure exponential factor
    double kinf_C_ = 0.0;
    /// @brief Troe F_c parameter
    double Fc_ = 0.6;
    /// @brief Troe N parameter
    double N_ = 1.0;
  };
}  // namespace micm
