// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file parameterized_function.hpp
/// @brief Device-safe POD for parameterizing a species concentration as a
///        linear combination of the standard Conditions fields.
///
/// Species::parameterize_ was previously a std::function, which cannot be
/// copied to device memory.  All existing MICM uses of parameterize_ are
/// linear combinations of {1, temperature, pressure, air_density}; storing
/// the coefficients directly gives a trivially-copyable POD that is safe to
/// deep-copy to a CUDA/HIP device via Kokkos.
///
/// For a nonlinear parameterization, extend ParameterizedFunction with more
/// coefficient fields (e.g. cross terms) rather than falling back to
/// std::function; anything reachable from device code must be POD.

#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

namespace micm
{
  /// @brief Linear parameterization
  ///        result = c0_ + c_T_ * T + c_P_ * P + c_rho_ * air_density
  struct ParameterizedFunction
  {
    Real c0_ = 0.0;
    Real c_T_ = 0.0;
    Real c_P_ = 0.0;
    Real c_rho_ = 0.0;
    Bool has_value_ = false;

    MICM_DEVICE_FUNCTION Real operator()(const Conditions& c) const
    {
      return c0_ + c_T_ * c.temperature_ + c_P_ * c.pressure_ + c_rho_ * c.air_density_;
    }

    MICM_DEVICE_FUNCTION Bool has_value() const
    {
      return has_value_;
    }
  };
}  // namespace micm
