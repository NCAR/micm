// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file reaction_rate_store.hpp
/// @brief Structure-of-arrays store for all reaction rate constant parameters.
///
/// Reactions must be stable-sorted by RateConstantTypeOrder before BuildFrom is
/// called (the SolverBuilder does this).  With sorted input each type occupies a
/// contiguous block and cumulative-size offsets make it easy to write batch
/// calculations directly into the correct locations in state.rate_constants_.  The
/// offsets are computed using the inline helper functions below.
///
/// **Two-phase calculation per step:**
///   Phase 1 — EvaluateCpuRates:  evaluate all LambdaRateConstantParameters entries
///                                  per grid cell, write results to state.rate_constants_.
///   Phase 2 — CalculateRates:    call batch analytic functions at the correct offsets,
///                                  then apply parameterized-species multipliers.

#define _USE_MATH_DEFINES

#include <micm/process/process.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/property_keys.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <math.h>
#include <type_traits>
#include <vector>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

namespace micm
{

  /// @brief Structure-of-arrays container for all reaction rate constant parameters.
  ///
  /// Requires that processes are sorted by RateConstantTypeOrder before BuildFrom is called.
  /// Each type occupies a contiguous block within state.rate_constants_[cell]; the block
  /// start is given by the inline offset helpers below.
  struct ReactionRateStore
  {
    // ----------------------------------------------------------------
    // Analytic types (GPU-safe parameter structs)
    // ----------------------------------------------------------------
    std::vector<ArrheniusRateConstantParameters> arrhenius;
    std::vector<TroeRateConstantParameters> troe;
    std::vector<TernaryChemicalActivationRateConstantParameters> ternary;
    std::vector<BranchedRateConstantParameters> branched;
    std::vector<TunnelingRateConstantParameters> tunneling;
    std::vector<TaylorSeriesRateConstantParameters> taylor;
    std::vector<ReversibleRateConstantParameters> reversible;

    // ----------------------------------------------------------------
    // Types using GPU-safe companion data structs
    // ----------------------------------------------------------------
    std::vector<UserDefinedRateConstantData> user_defined;
    std::vector<SurfaceRateConstantData> surface;

    // ----------------------------------------------------------------
    // Contiguous-block offsets into state.rate_constants_[cell]
    // ----------------------------------------------------------------
    std::size_t troe_offset() const
    {
      return arrhenius.size();
    }
    std::size_t ternary_offset() const
    {
      return troe_offset() + troe.size();
    }
    std::size_t branched_offset() const
    {
      return ternary_offset() + ternary.size();
    }
    std::size_t tunneling_offset() const
    {
      return branched_offset() + branched.size();
    }
    std::size_t taylor_offset() const
    {
      return tunneling_offset() + tunneling.size();
    }
    std::size_t reversible_offset() const
    {
      return taylor_offset() + taylor.size();
    }
    std::size_t user_defined_offset() const
    {
      return reversible_offset() + reversible.size();
    }
    std::size_t surface_offset() const
    {
      return user_defined_offset() + user_defined.size();
    }
    std::size_t lambda_offset() const
    {
      return surface_offset() + surface.size();
    }

    // ----------------------------------------------------------------
    // CPU-only lambda entries
    // ----------------------------------------------------------------

    /// @brief One entry per LambdaRateConstantParameters reaction.
    ///        Stores a non-owning pointer into the Solver's processes_ so that
    ///        mutations via Solver::GetLambdaRateConstantByName propagate immediately.
    struct LambdaEntry
    {
      /// @brief Non-owning pointer; valid for the lifetime of the owning Solver.
      LambdaRateConstantParameters* source;
      /// @brief Column index in state.rate_constants_[cell] for this reaction.
      std::size_t rc_index;
    };
    std::vector<LambdaEntry> lambda_entries;

    // ----------------------------------------------------------------
    // Parameterized-species multipliers
    // ----------------------------------------------------------------

    /// @brief Applied as a final pass after all batch calculations.
    ///        One entry per reaction that has at least one parameterized reactant.
    struct ParameterizedMultiplier
    {
      std::function<double(const Conditions&)> evaluate;
      std::size_t rc_index;
    };
    std::vector<ParameterizedMultiplier> parameterized_multipliers;

    // ================================================================
    // Factory
    // ================================================================

    /// @brief Build a ReactionRateStore from a sorted process list.
    ///
    ///        Requires that processes are stable-sorted by RateConstantTypeOrder
    ///        (done by SolverBuilder) before this call so the type blocks are contiguous
    ///        and the offset helpers return correct values.
    ///
    /// @param processes  Non-const ref to the solver's process list (allows non-const
    ///                   pointers into LambdaRateConstantParameters for runtime mutation)
    /// @return           Populated store ready for two-phase calculation
    static ReactionRateStore BuildFrom(std::vector<Process>& processes)
    {
      ReactionRateStore store;
      std::size_t rc_index = 0;
      std::size_t custom_param_off = 0;

      for (auto& process : processes)
      {
        ChemicalReaction* reaction = std::get_if<ChemicalReaction>(&process.process_);
        if (!reaction)
          continue;

        // Collect parameterized-reactant multiplier for this reaction.
        {
          std::vector<std::function<double(const Conditions&)>> param_funcs;
          for (const auto& reactant : reaction->reactants_)
            if (reactant.IsParameterized())
              param_funcs.push_back(reactant.parameterize_);

          if (!param_funcs.empty())
          {
            store.parameterized_multipliers.push_back(
                { [pf = std::move(param_funcs)](const Conditions& cond)
                  {
                    double val = 1.0;
                    for (const auto& f : pf)
                      val *= f(cond);
                    return val;
                  },
                  rc_index });
          }
        }

        std::size_t n_custom = 0;

        RateConstantVariant& rc = reaction->rate_constant_;

        if (auto* p = std::get_if<ArrheniusRateConstantParameters>(&rc))
        {
          store.arrhenius.push_back(*p);
        }
        else if (auto* p = std::get_if<TroeRateConstantParameters>(&rc))
        {
          store.troe.push_back(*p);
        }
        else if (auto* p = std::get_if<TernaryChemicalActivationRateConstantParameters>(&rc))
        {
          store.ternary.push_back(*p);
        }
        else if (auto* p = std::get_if<BranchedRateConstantParameters>(&rc))
        {
          BranchedRateConstantParameters params = *p;
          // Compute precomputed derived fields (same formulae as the old BranchedRateConstant ctor)
          params.k0_ = 2.0e-22 * constants::AVOGADRO_CONSTANT * 1.0e-6 * std::exp(static_cast<double>(params.n_));
          double air_ref = 2.45e19 / constants::AVOGADRO_CONSTANT * 1.0e6;
          double a = params.k0_ * air_ref;
          double b = 0.43 * std::pow(293.0 / 298.0, -8.0);
          double A_val = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
          params.z_ = A_val * (1.0 - params.a0_) / params.a0_;
          store.branched.push_back(params);
        }
        else if (auto* p = std::get_if<TunnelingRateConstantParameters>(&rc))
        {
          store.tunneling.push_back(*p);
        }
        else if (auto* p = std::get_if<TaylorSeriesRateConstantParameters>(&rc))
        {
          store.taylor.push_back(*p);
        }
        else if (auto* p = std::get_if<ReversibleRateConstantParameters>(&rc))
        {
          store.reversible.push_back(*p);
        }
        else if (auto* p = std::get_if<UserDefinedRateConstantParameters>(&rc))
        {
          UserDefinedRateConstantData data;
          data.scaling_factor_ = p->scaling_factor_;
          data.custom_param_index_ = custom_param_off;
          store.user_defined.push_back(data);
          n_custom = 1;
        }
        else if (auto* p = std::get_if<SurfaceRateConstantParameters>(&rc))
        {
          SurfaceRateConstantData data;
          if (!p->phase_species_.diffusion_coefficient_.has_value())
            throw MicmException(
                MicmSeverity::Error,
                MICM_ERROR_CATEGORY_SPECIES,
                MICM_SPECIES_ERROR_CODE_PROPERTY_NOT_FOUND,
                "Diffusion coefficient for species '" + p->phase_species_.species_.name_ + "' is not defined");
          data.diffusion_coefficient_ = p->phase_species_.diffusion_coefficient_.value();
          double mw = p->phase_species_.species_.GetProperty<double>(property_keys::MOLECULAR_WEIGHT);
          data.mean_free_speed_factor_ = 8.0 * constants::GAS_CONSTANT / (M_PI * mw);
          data.reaction_probability_ = p->reaction_probability_;
          data.custom_param_base_index_ = custom_param_off;
          store.surface.push_back(data);
          n_custom = 2;
        }
        else if (auto* p = std::get_if<LambdaRateConstantParameters>(&rc))
        {
          store.lambda_entries.push_back({ p, rc_index });
        }

        custom_param_off += n_custom;
        ++rc_index;
      }

      return store;
    }

    // ================================================================
    // Phase 1 — CPU-only lambda evaluation
    // ================================================================

    /// @brief Evaluate all lambda rate constants and write results to
    ///        state.rate_constants_.  Must be called every step before CalculateRates.
    template<class StatePolicy>
    static void EvaluateCpuRates(const ReactionRateStore& store, StatePolicy& state)
    {
      if (store.lambda_entries.empty())
        return;

      using DenseMatrixPolicy = typename StatePolicy::DenseMatrixPolicyType;

      if constexpr (VectorizableDense<DenseMatrixPolicy>)
      {
        auto& v_rc = state.rate_constants_.AsVector();
        constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();

        for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
        {
          const std::size_t rc_base = i_group * state.rate_constants_.GroupSize();
          const std::size_t n_local = std::min(L, state.rate_constants_.NumRows() - i_group * L);

          for (std::size_t i_cell = 0; i_cell < n_local; ++i_cell)
          {
            const auto& cond = state.conditions_[i_group * L + i_cell];
            for (const auto& entry : store.lambda_entries)
              v_rc[rc_base + entry.rc_index * L + i_cell] = entry.source->lambda_function_(cond);
          }
        }
      }
      else
      {
        for (std::size_t i_cell = 0; i_cell < state.rate_constants_.NumRows(); ++i_cell)
        {
          const auto& cond = state.conditions_[i_cell];
          for (const auto& entry : store.lambda_entries)
            state.rate_constants_[i_cell][entry.rc_index] = entry.source->lambda_function_(cond);
        }
      }
    }

    // ================================================================
    // Phase 2 — batch analytic calculation
    // ================================================================

    /// @brief Calculate all analytic rate constants and write them into
    ///        state.rate_constants_ at the correct type-block offsets.
    ///        Lambda entries (written by EvaluateCpuRates) are untouched.
    ///        Parameterized-species multipliers are applied last.
    template<class StatePolicy>
    static void CalculateRates(const ReactionRateStore& store, StatePolicy& state)
    {
      using DenseMatrixPolicy = typename StatePolicy::DenseMatrixPolicyType;
      if constexpr (VectorizableDense<DenseMatrixPolicy>)
        CalculateRatesVectorized<DenseMatrixPolicy>(store, state);
      else
        CalculateRatesScalar(store, state);
    }

   private:
    // ----------------------------------------------------------------
    // Scalar implementation
    // ----------------------------------------------------------------

    template<class StatePolicy>
    static void CalculateRatesScalar(const ReactionRateStore& store, StatePolicy& state)
    {
      const std::size_t n_cells = state.rate_constants_.NumRows();

      for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
      {
        const auto& cond = state.conditions_[i_cell];
        const std::vector<double> cp_row = state.custom_rate_parameters_[i_cell];
        const double* cp = cp_row.data();

        // Helper: call batch function, write directly at offset+k
        auto run = [&](auto batch_fn, std::size_t n, std::size_t offset)
        {
          if (n == 0)
            return;
          std::vector<double> buf(n);
          batch_fn(buf.data());
          for (std::size_t k = 0; k < n; ++k)
            state.rate_constants_[i_cell][offset + k] = buf[k];
        };

        run([&](double* buf)
            { CalculateArrhenius(store.arrhenius.data(), store.arrhenius.size(), cond.temperature_, cond.pressure_, buf); },
            store.arrhenius.size(),
            0);

        run([&](double* buf)
            { CalculateTroe(store.troe.data(), store.troe.size(), cond.temperature_, cond.air_density_, buf); },
            store.troe.size(),
            store.troe_offset());

        run(
            [&](double* buf)
            {
              CalculateTernaryChemicalActivation(
                  store.ternary.data(), store.ternary.size(), cond.temperature_, cond.air_density_, buf);
            },
            store.ternary.size(),
            store.ternary_offset());

        run([&](double* buf)
            { CalculateBranched(store.branched.data(), store.branched.size(), cond.temperature_, cond.air_density_, buf); },
            store.branched.size(),
            store.branched_offset());

        run([&](double* buf) { CalculateTunneling(store.tunneling.data(), store.tunneling.size(), cond.temperature_, buf); },
            store.tunneling.size(),
            store.tunneling_offset());

        run([&](double* buf)
            { CalculateTaylorSeries(store.taylor.data(), store.taylor.size(), cond.temperature_, cond.pressure_, buf); },
            store.taylor.size(),
            store.taylor_offset());

        run([&](double* buf)
            { CalculateReversible(store.reversible.data(), store.reversible.size(), cond.temperature_, buf); },
            store.reversible.size(),
            store.reversible_offset());

        run([&](double* buf) { CalculateUserDefined(store.user_defined.data(), store.user_defined.size(), cp, buf); },
            store.user_defined.size(),
            store.user_defined_offset());

        run([&](double* buf) { CalculateSurface(store.surface.data(), store.surface.size(), cond.temperature_, cp, buf); },
            store.surface.size(),
            store.surface_offset());

        for (const auto& mult : store.parameterized_multipliers)
          state.rate_constants_[i_cell][mult.rc_index] *= mult.evaluate(cond);
      }
    }

    // ----------------------------------------------------------------
    // Vectorizable implementation
    // ----------------------------------------------------------------

    template<class DenseMatrixPolicy, class StatePolicy>
    static void CalculateRatesVectorized(const ReactionRateStore& store, StatePolicy& state)
    {
      auto& v_rc = state.rate_constants_.AsVector();
      const auto& v_cp = state.custom_rate_parameters_.AsVector();
      constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
      const std::size_t n_cp = state.custom_rate_parameters_.NumColumns();

      for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
      {
        const std::size_t rc_base = i_group * state.rate_constants_.GroupSize();
        const std::size_t cp_base = i_group * state.custom_rate_parameters_.GroupSize();
        const std::size_t n_local = std::min(L, state.rate_constants_.NumRows() - i_group * L);

        for (std::size_t i_cell = 0; i_cell < n_local; ++i_cell)
        {
          const auto& cond = state.conditions_[i_group * L + i_cell];

          // Gather per-cell custom params from interleaved storage
          std::vector<double> cp(n_cp);
          for (std::size_t j = 0; j < n_cp; ++j)
            cp[j] = v_cp[cp_base + j * L + i_cell];

          // Helper: call batch function, write at offset+k in interleaved layout
          auto run = [&](auto batch_fn, std::size_t n, std::size_t offset)
          {
            if (n == 0)
              return;
            std::vector<double> buf(n);
            batch_fn(buf.data());
            for (std::size_t k = 0; k < n; ++k)
              v_rc[rc_base + (offset + k) * L + i_cell] = buf[k];
          };

          run(
              [&](double* buf)
              {
                CalculateArrhenius(store.arrhenius.data(), store.arrhenius.size(), cond.temperature_, cond.pressure_, buf);
              },
              store.arrhenius.size(),
              0);

          run([&](double* buf)
              { CalculateTroe(store.troe.data(), store.troe.size(), cond.temperature_, cond.air_density_, buf); },
              store.troe.size(),
              store.troe_offset());

          run(
              [&](double* buf)
              {
                CalculateTernaryChemicalActivation(
                    store.ternary.data(), store.ternary.size(), cond.temperature_, cond.air_density_, buf);
              },
              store.ternary.size(),
              store.ternary_offset());

          run(
              [&](double* buf)
              {
                CalculateBranched(store.branched.data(), store.branched.size(), cond.temperature_, cond.air_density_, buf);
              },
              store.branched.size(),
              store.branched_offset());

          run([&](double* buf)
              { CalculateTunneling(store.tunneling.data(), store.tunneling.size(), cond.temperature_, buf); },
              store.tunneling.size(),
              store.tunneling_offset());

          run([&](double* buf)
              { CalculateTaylorSeries(store.taylor.data(), store.taylor.size(), cond.temperature_, cond.pressure_, buf); },
              store.taylor.size(),
              store.taylor_offset());

          run([&](double* buf)
              { CalculateReversible(store.reversible.data(), store.reversible.size(), cond.temperature_, buf); },
              store.reversible.size(),
              store.reversible_offset());

          run([&](double* buf)
              { CalculateUserDefined(store.user_defined.data(), store.user_defined.size(), cp.data(), buf); },
              store.user_defined.size(),
              store.user_defined_offset());

          run([&](double* buf)
              { CalculateSurface(store.surface.data(), store.surface.size(), cond.temperature_, cp.data(), buf); },
              store.surface.size(),
              store.surface_offset());

          for (const auto& mult : store.parameterized_multipliers)
            v_rc[rc_base + mult.rc_index * L + i_cell] *= mult.evaluate(cond);
        }
      }
    }
  };

}  // namespace micm
