// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file rate_constant_functions.hpp
/// @brief Calculation functions for all supported rate constant types.
///
/// Each function processes a contiguous block of reactions of a single type for
/// one grid cell.  The loop is inside the function so the compiler can
/// SIMD-vectorize the inner pass; on GPU the same functions run per-thread via
/// constexpr (callable from device code with --expt-relaxed-constexpr).
///
/// **Output pointer contract:**
/// The caller must position `output` at the correct type-group offset within
/// `rate_constants[cell]`.  ReactionRateConstantStore provides inline offset helpers
/// (e.g. troe_offset(), ternary_offset()) for this purpose.
///
/// **CPU/GPU split:**
/// Analytic functions (CalculateArrhenius, CalculateTroe, …) are safe on both
/// CPU and GPU.  CalculateUserDefined and CalculateSurface read from
/// `custom_params`, which must be the device pointer when called from a CUDA kernel.

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
#  define MICM_CONSTEXPR constexpr
#else
#  define MICM_CONSTEXPR
#endif

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace micm
{

  /// @brief Calculate Arrhenius rate constants for a contiguous block of reactions.
  ///        k = A * exp(C/T) * (T/D)^B * (1 + E*P)
  /// @param params      Array of Arrhenius parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param pressure    Pressure [Pa]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateArrhenius(
      const ArrheniusRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double pressure,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = params[i].A_ * std::exp(params[i].C_ / temperature) *
                  std::pow(temperature / params[i].D_, params[i].B_) *
                  (1.0 + params[i].E_ * pressure);
  }

  /// @brief Shared falloff kernel for Troe and TernaryChemicalActivation.
  ///        result = k0 * numerator_scale / (1 + ratio) * Fc^(N/(N + log10(ratio)^2))
  ///        Troe passes air_density as numerator_scale; Ternary passes 1.0.
  template<class FalloffParams>
  MICM_CONSTEXPR inline double FalloffKernel(
      const FalloffParams& p,
      double temperature,
      double air_density,
      double numerator_scale)
  {
    double k0    = p.k0_A_ * std::exp(p.k0_C_ / temperature) * std::pow(temperature / 300.0, p.k0_B_);
    double kinf  = p.kinf_A_ * std::exp(p.kinf_C_ / temperature) * std::pow(temperature / 300.0, p.kinf_B_);
    double ratio = k0 * air_density / kinf;
    return k0 * numerator_scale / (1.0 + ratio) *
           std::pow(p.Fc_, p.N_ / (p.N_ + std::pow(std::log10(ratio), 2.0)));
  }

  /// @brief Calculate Troe rate constants for a contiguous block of reactions.
  /// @param params      Array of Troe parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param air_density Air number density [mol m-3]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateTroe(
      const TroeRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double air_density,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = FalloffKernel(params[i], temperature, air_density, air_density);
  }

  /// @brief Calculate Ternary Chemical Activation rate constants for a contiguous block of reactions.
  /// @param params      Array of TernaryChemicalActivation parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param air_density Air number density [mol m-3]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateTernaryChemicalActivation(
      const TernaryChemicalActivationRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double air_density,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = FalloffKernel(params[i], temperature, air_density, 1.0);
  }

  /// @brief Calculate Tunneling rate constants for a contiguous block of reactions.
  ///        k = A * exp(-B/T + C/T^3)
  /// @param params      Array of Tunneling parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateTunneling(
      const TunnelingRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = params[i].A_ *
                  std::exp(-params[i].B_ / temperature +
                           params[i].C_ / (temperature * temperature * temperature));
  }

  /// @brief Calculate Branched rate constants for a contiguous block of reactions.
  ///        Requires params[i].k0_ and params[i].z_ to be precomputed by
  ///        ReactionRateConstantStore::BuildFrom before calling this function.
  /// @param params      Array of Branched parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param air_density Air number density [mol m-3]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateBranched(
      const BranchedRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double air_density,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      double a     = params[i].k0_ * air_density;
      double b     = 0.43 * std::pow(temperature / 298.0, -8.0);
      double A_val = a / (1.0 + a / b) *
                     std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
      double pre = params[i].X_ * std::exp(-params[i].Y_ / temperature);
      output[i]  = (params[i].branch_ == BranchedRateConstantParameters::Branch::Alkoxy)
                       ? pre * (params[i].z_ / (params[i].z_ + A_val))
                       : pre * (A_val / (A_val + params[i].z_));
    }
  }

  /// @brief Calculate Taylor Series rate constants for a contiguous block of reactions.
  ///        k = (sum_{j=0}^{n-1} c_j * T^j) * A * exp(C/T) * (T/D)^B * (1 + E*P)
  /// @param params      Array of TaylorSeries parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param pressure    Pressure [Pa]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateTaylorSeries(
      const TaylorSeriesRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double pressure,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
    {
      double poly  = 0.0;
      double t_pow = 1.0;
      for (std::size_t j = 0; j < params[i].n_coefficients_; ++j, t_pow *= temperature)
        poly += params[i].coefficients_[j] * t_pow;
      output[i] = poly * params[i].A_ * std::exp(params[i].C_ / temperature) *
                  std::pow(temperature / params[i].D_, params[i].B_) *
                  (1.0 + params[i].E_ * pressure);
    }
  }

  /// @brief Calculate Reversible rate constants for a contiguous block of reactions.
  ///        k = A * exp(C/T) * k_r
  /// @param params      Array of Reversible parameters, length n
  /// @param n           Number of reactions
  /// @param temperature Temperature [K]
  /// @param output      Destination array of length n
  MICM_CONSTEXPR inline void CalculateReversible(
      const ReversibleRateConstantParameters* params,
      std::size_t n,
      double temperature,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = params[i].A_ * std::exp(params[i].C_ / temperature) * params[i].k_r_;
  }

  /// @brief Apply user-defined rate constants from per-grid-cell custom parameters.
  ///        output[i] = custom_params[params[i].custom_param_index_] * params[i].scaling_factor_
  /// @param params        Array of UserDefined data, length n
  /// @param n             Number of reactions
  /// @param custom_params Row of custom_rate_parameters_ for this grid cell
  /// @param output        Destination array of length n
  MICM_CONSTEXPR inline void CalculateUserDefined(
      const UserDefinedRateConstantData* params,
      std::size_t n,
      const double* custom_params,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = custom_params[params[i].custom_param_index_] * params[i].scaling_factor_;
  }

  /// @brief Calculate one surface rate constant given pre-fetched aerosol parameters.
  /// @param radius    Aerosol effective radius [m]
  /// @param num_conc  Particle number concentration [# m-3]
  MICM_CONSTEXPR inline double CalculateSurfaceOne(
      const SurfaceRateConstantData& p,
      double temperature,
      double radius,
      double num_conc)
  {
    double mean_free_speed = std::sqrt(p.mean_free_speed_factor_ * temperature);
    return 4.0 * num_conc * M_PI * radius * radius /
           (radius / p.diffusion_coefficient_ + 4.0 / (mean_free_speed * p.reaction_probability_));
  }

  /// @brief Calculate surface rate constants for n reactions.
  ///        Reads radius and num_conc from custom_params[custom_param_base_index_] and +1.
  MICM_CONSTEXPR inline void CalculateSurface(
      const SurfaceRateConstantData* params,
      std::size_t n,
      double temperature,
      const double* custom_params,
      double* output)
  {
    for (std::size_t i = 0; i < n; ++i)
      output[i] = CalculateSurfaceOne(
          params[i],
          temperature,
          custom_params[params[i].custom_param_base_index_],
          custom_params[params[i].custom_param_base_index_ + 1]);
  }

}  // namespace micm
