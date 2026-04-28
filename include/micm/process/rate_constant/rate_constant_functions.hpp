// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file rate_constant_functions.hpp
/// @brief Calculation functions for all supported rate constant types.
///
/// Each function takes a single parameter struct and returns a single double.
/// CPU: called inside ForEachRow loops (one reaction at a time, all cells).
/// GPU: called per-thread inside CalculateRatesForThread (one cell per thread).
///
/// All functions are MICM_CONSTEXPR so they are callable from device code
/// with --expt-relaxed-constexpr.

#define _USE_MATH_DEFINES
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/taylor_series_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>

#include <cmath>
#include <cstddef>
#include <math.h>

// constexpr <cmath> functions (std::exp, std::pow, std::sqrt, etc.) require C++23
// and compiler support (P1383R2). Guard so MSVC and other compilers that haven't
// implemented this yet still produce valid code.
// __CUDACC__: nvcc always needs constexpr so --expt-relaxed-constexpr can call
// these from device code, even when the host compiler lacks P1383R2.
#if defined(__cpp_lib_constexpr_cmath) || defined(__CUDACC__)
  #define MICM_CONSTEXPR constexpr
#else
  #define MICM_CONSTEXPR
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

namespace micm
{

  /// @brief Shared falloff kernel for Troe and TernaryChemicalActivation.
  ///        result = k0 * numerator_scale / (1 + ratio) * Fc^(N/(N + log10(ratio)^2))
  ///        Troe passes air_density as numerator_scale; Ternary passes 1.0.
  template<class FalloffParams>
  MICM_CONSTEXPR inline double
  FalloffKernel(const FalloffParams& p, double temperature, double air_density, double numerator_scale)
  {
    double k0 = p.k0_A_ * std::exp(p.k0_C_ / temperature) * std::pow(temperature / 300.0, p.k0_B_);
    double kinf = p.kinf_A_ * std::exp(p.kinf_C_ / temperature) * std::pow(temperature / 300.0, p.kinf_B_);
    double ratio = k0 * air_density / kinf;
    return k0 * numerator_scale / (1.0 + ratio) * std::pow(p.Fc_, p.N_ / (p.N_ + std::pow(std::log10(ratio), 2.0)));
  }

  /// @brief Calculate Arrhenius rate constant.
  ///        k = A * exp(C/T) * (T/D)^B * (1 + E*P)
  MICM_CONSTEXPR inline double
  CalculateArrhenius(const ArrheniusRateConstantParameters& p, double temperature, double pressure)
  {
    return p.A_ * std::exp(p.C_ / temperature) * std::pow(temperature / p.D_, p.B_) * (1.0 + p.E_ * pressure);
  }

  /// @brief Calculate Troe rate constant.
  MICM_CONSTEXPR inline double CalculateTroe(const TroeRateConstantParameters& p, double temperature, double air_density)
  {
    return FalloffKernel(p, temperature, air_density, air_density);
  }

  /// @brief Calculate Ternary Chemical Activation rate constant.
  MICM_CONSTEXPR inline double CalculateTernaryChemicalActivation(
      const TernaryChemicalActivationRateConstantParameters& p,
      double temperature,
      double air_density)
  {
    return FalloffKernel(p, temperature, air_density, 1.0);
  }

  /// @brief Calculate Tunneling rate constant.
  ///        k = A * exp(-B/T + C/T^3)
  MICM_CONSTEXPR inline double CalculateTunneling(const TunnelingRateConstantParameters& p, double temperature)
  {
    return p.A_ * std::exp(-p.B_ / temperature + p.C_ / (temperature * temperature * temperature));
  }

  /// @brief Calculate Branched rate constant.
  ///        Requires p.k0_ and p.z_ to be precomputed by ReactionRateConstantStore::BuildFrom.
  MICM_CONSTEXPR inline double
  CalculateBranched(const BranchedRateConstantParameters& p, double temperature, double air_density)
  {
    double a = p.k0_ * air_density;
    double b = 0.43 * std::pow(temperature / 298.0, -8.0);
    double A_val = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
    double pre = p.X_ * std::exp(-p.Y_ / temperature);
    return (p.branch_ == BranchedRateConstantParameters::Branch::Alkoxy) ? pre * (p.z_ / (p.z_ + A_val))
                                                                         : pre * (A_val / (A_val + p.z_));
  }

  /// @brief Calculate Taylor Series rate constant.
  ///        k = (sum_{j=0}^{n-1} c_j * T^j) * A * exp(C/T) * (T/D)^B * (1 + E*P)
  MICM_CONSTEXPR inline double
  CalculateTaylorSeries(const TaylorSeriesRateConstantParameters& p, double temperature, double pressure)
  {
    double poly = 0.0;
    double t_pow = 1.0;
    for (std::size_t j = 0; j < p.n_coefficients_; ++j, t_pow *= temperature)
      poly += p.coefficients_[j] * t_pow;
    return poly * p.A_ * std::exp(p.C_ / temperature) * std::pow(temperature / p.D_, p.B_) * (1.0 + p.E_ * pressure);
  }

  /// @brief Calculate Reversible rate constant.
  ///        k = A * exp(C/T) * k_r
  MICM_CONSTEXPR inline double CalculateReversible(const ReversibleRateConstantParameters& p, double temperature)
  {
    return p.A_ * std::exp(p.C_ / temperature) * p.k_r_;
  }

  /// @brief Calculate user-defined rate constant.
  ///        k = custom_param_value * scaling_factor
  MICM_CONSTEXPR inline double CalculateUserDefined(const UserDefinedRateConstantData& p, double custom_param_value)
  {
    return custom_param_value * p.scaling_factor_;
  }

  /// @brief Calculate one surface rate constant given pre-fetched aerosol parameters.
  /// @param radius    Aerosol effective radius [m]
  /// @param num_conc  Particle number concentration [# m-3]
  MICM_CONSTEXPR inline double
  CalculateSurfaceOne(const SurfaceRateConstantData& p, double temperature, double radius, double num_conc)
  {
    double mean_free_speed = std::sqrt(p.mean_free_speed_factor_ * temperature);
    return 4.0 * num_conc * M_PI * radius * radius /
           (radius / p.diffusion_coefficient_ + 4.0 / (mean_free_speed * p.reaction_probability_));
  }

}  // namespace micm
