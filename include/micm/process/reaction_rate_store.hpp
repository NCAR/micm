// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file reaction_rate_store.hpp
/// @brief Structure-of-arrays store for all reaction rate constant parameters.
///
/// Reactions must be stable-sorted by RateConstantTypeOrder before BuildFrom is called
/// (the SolverBuilder does this).  Each type occupies a contiguous block; offset helpers
/// below give the start of each block within state.rate_constants_[cell].
///
/// Each step: call EvaluateCpuRateConstants (lambda entries), then CpuCalculateRateConstants (analytic types).

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

  /// @brief Structure-of-arrays store for all reaction rate constant parameters.
  ///
  /// Processes must be sorted by RateConstantTypeOrder before BuildFrom is called.
  struct ReactionRateConstantStore
  {
    // ----------------------------------------------------------------
    // Analytic types (GPU-safe parameter structs)
    // ----------------------------------------------------------------
    std::vector<ArrheniusRateConstantParameters> arrhenius_;
    std::vector<TroeRateConstantParameters> troe_;
    std::vector<TernaryChemicalActivationRateConstantParameters> ternary_;
    std::vector<BranchedRateConstantParameters> branched_;
    std::vector<TunnelingRateConstantParameters> tunneling_;
    std::vector<TaylorSeriesRateConstantParameters> taylor_;
    std::vector<ReversibleRateConstantParameters> reversible_;

    // ----------------------------------------------------------------
    // Types using GPU-safe companion data structs
    // ----------------------------------------------------------------
    std::vector<UserDefinedRateConstantData> user_defined_;
    std::vector<SurfaceRateConstantData> surface_;

    // ----------------------------------------------------------------
    // Contiguous-block offsets into state.rate_constants_[cell]
    // ----------------------------------------------------------------
    std::size_t troe_offset() const
    {
      return arrhenius_.size();
    }
    std::size_t ternary_offset() const
    {
      return troe_offset() + troe_.size();
    }
    std::size_t branched_offset() const
    {
      return ternary_offset() + ternary_.size();
    }
    std::size_t tunneling_offset() const
    {
      return branched_offset() + branched_.size();
    }
    std::size_t taylor_offset() const
    {
      return tunneling_offset() + tunneling_.size();
    }
    std::size_t reversible_offset() const
    {
      return taylor_offset() + taylor_.size();
    }
    std::size_t user_defined_offset() const
    {
      return reversible_offset() + reversible_.size();
    }
    std::size_t surface_offset() const
    {
      return user_defined_offset() + user_defined_.size();
    }
    std::size_t lambda_offset() const
    {
      return surface_offset() + surface_.size();
    }

    // ----------------------------------------------------------------
    // CPU-only lambda entries
    // ----------------------------------------------------------------

    struct LambdaEntry
    {
      LambdaRateConstantParameters* source;  ///< Non-owning; valid for the lifetime of the owning Solver.
      std::size_t rc_index;                  ///< Column index in state.rate_constants_[cell].
    };
    std::vector<LambdaEntry> lambda_entries_;

    // ----------------------------------------------------------------
    // Parameterized-species multipliers
    // ----------------------------------------------------------------

    /// @brief One entry per reaction with at least one parameterized reactant.
    struct ParameterizedMultiplier
    {
      std::function<double(const Conditions&)> evaluate;
      std::size_t rc_index;
    };
    std::vector<ParameterizedMultiplier> parameterized_multipliers_;

    // ================================================================
    // Factory
    // ================================================================

    /// @brief Build a ReactionRateConstantStore from a sorted process list.
    /// @param processes  Non-const ref so LambdaRateConstantParameters pointers remain mutable at runtime.
    static ReactionRateConstantStore BuildFrom(std::vector<Process>& processes)
    {
      ReactionRateConstantStore store;
      std::size_t rc_index = 0;
      std::size_t custom_param_off = 0;

      for (auto& process : processes)
      {
        ChemicalReaction* reaction = std::get_if<ChemicalReaction>(&process.process_);
        if (!reaction)
          continue;

        {  // parameterized-reactant multiplier
          std::vector<std::function<double(const Conditions&)>> param_funcs;
          for (const auto& reactant : reaction->reactants_)
            if (reactant.IsParameterized())
              param_funcs.push_back(reactant.parameterize_);

          if (!param_funcs.empty())
          {
            store.parameterized_multipliers_.push_back(
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
          store.arrhenius_.push_back(*p);
        }
        else if (auto* p = std::get_if<TroeRateConstantParameters>(&rc))
        {
          store.troe_.push_back(*p);
        }
        else if (auto* p = std::get_if<TernaryChemicalActivationRateConstantParameters>(&rc))
        {
          store.ternary_.push_back(*p);
        }
        else if (auto* p = std::get_if<BranchedRateConstantParameters>(&rc))
        {
          BranchedRateConstantParameters params = *p;
          // Pre-compute derived fields needed by CalculateBranched
          params.k0_ = 2.0e-22 * constants::AVOGADRO_CONSTANT * 1.0e-6 * std::exp(static_cast<double>(params.n_));
          double air_ref = 2.45e19 / constants::AVOGADRO_CONSTANT * 1.0e6;
          double a = params.k0_ * air_ref;
          double b = 0.43 * std::pow(293.0 / 298.0, -8.0);
          double A_val = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
          params.z_ = A_val * (1.0 - params.a0_) / params.a0_;
          store.branched_.push_back(params);
        }
        else if (auto* p = std::get_if<TunnelingRateConstantParameters>(&rc))
        {
          store.tunneling_.push_back(*p);
        }
        else if (auto* p = std::get_if<TaylorSeriesRateConstantParameters>(&rc))
        {
          store.taylor_.push_back(*p);
        }
        else if (auto* p = std::get_if<ReversibleRateConstantParameters>(&rc))
        {
          store.reversible_.push_back(*p);
        }
        else if (auto* p = std::get_if<UserDefinedRateConstantParameters>(&rc))
        {
          UserDefinedRateConstantData data;
          data.scaling_factor_ = p->scaling_factor_;
          data.custom_param_index_ = custom_param_off;
          store.user_defined_.push_back(data);
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
          store.surface_.push_back(data);
          n_custom = 2;
        }
        else if (auto* p = std::get_if<LambdaRateConstantParameters>(&rc))
        {
          store.lambda_entries_.push_back({ p, rc_index });
        }

        custom_param_off += n_custom;
        ++rc_index;
      }

      return store;
    }

    /// @brief Evaluate all lambda rate constants into state.rate_constants_.
    ///        Must be called each step before CpuCalculateRateConstants.
    template<class StatePolicy>
    static void EvaluateCpuRateConstants(const ReactionRateConstantStore& store, StatePolicy& state)
    {
      if (store.lambda_entries_.empty())
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
            for (const auto& entry : store.lambda_entries_)
              v_rc[rc_base + entry.rc_index * L + i_cell] = entry.source->lambda_function_(cond);
          }
        }
      }
      else
      {
        for (std::size_t i_cell = 0; i_cell < state.rate_constants_.NumRows(); ++i_cell)
        {
          const auto& cond = state.conditions_[i_cell];
          for (const auto& entry : store.lambda_entries_)
            state.rate_constants_[i_cell][entry.rc_index] = entry.source->lambda_function_(cond);
        }
      }
    }

    /// @brief Calculate all analytic rate constants into state.rate_constants_.
    ///        Lambda entries are untouched; parameterized multipliers applied last.
    template<class StatePolicy>
    static void CpuCalculateRateConstants(const ReactionRateConstantStore& store, StatePolicy& state)
    {
      using DenseMatrixPolicy = typename StatePolicy::DenseMatrixPolicyType;
      if constexpr (VectorizableDense<DenseMatrixPolicy>)
        CpuCalculateRateConstantsVectorized<DenseMatrixPolicy>(store, state);
      else
        CpuCalculateRateConstantsScalar(store, state);
    }

   private:
    /// @brief Calculate all analytic rate constants for one cell into a contiguous buffer.
    /// @param cp   Custom rate parameters for this cell (contiguous)
    /// @param out  Output buffer, size >= lambda_offset()
    static void CpuCalculateOneCellRateConstants(
        const ReactionRateConstantStore& store,
        double temperature,
        double pressure,
        double air_density,
        const double* cp,
        double* out)
    {
      CalculateArrhenius(store.arrhenius_.data(), store.arrhenius_.size(), temperature, pressure, out);
      CalculateTroe(store.troe_.data(), store.troe_.size(), temperature, air_density, out + store.troe_offset());
      CalculateTernaryChemicalActivation(
          store.ternary_.data(), store.ternary_.size(), temperature, air_density, out + store.ternary_offset());
      CalculateBranched(
          store.branched_.data(), store.branched_.size(), temperature, air_density, out + store.branched_offset());
      CalculateTunneling(store.tunneling_.data(), store.tunneling_.size(), temperature, out + store.tunneling_offset());
      CalculateTaylorSeries(
          store.taylor_.data(), store.taylor_.size(), temperature, pressure, out + store.taylor_offset());
      CalculateReversible(store.reversible_.data(), store.reversible_.size(), temperature, out + store.reversible_offset());
      CalculateUserDefined(store.user_defined_.data(), store.user_defined_.size(), cp, out + store.user_defined_offset());
      CalculateSurface(store.surface_.data(), store.surface_.size(), temperature, cp, out + store.surface_offset());
    }

    template<class StatePolicy>
    static void CpuCalculateRateConstantsScalar(const ReactionRateConstantStore& store, StatePolicy& state)
    {
      const std::size_t n_cells = state.rate_constants_.NumRows();
      const std::size_t n_rc = state.rate_constants_.NumColumns();
      const std::size_t n_cp = state.custom_rate_parameters_.NumColumns();
      double* rc_data = state.rate_constants_.AsVector().data();
      const double* cp_data = state.custom_rate_parameters_.AsVector().data();

      for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
      {
        const auto& cond = state.conditions_[i_cell];
        double* rc_row = rc_data + i_cell * n_rc;
        const double* cp_row = cp_data + i_cell * n_cp;

        CpuCalculateOneCellRateConstants(store, cond.temperature_, cond.pressure_, cond.air_density_, cp_row, rc_row);

        for (const auto& mult : store.parameterized_multipliers_)
          rc_row[mult.rc_index] *= mult.evaluate(cond);
      }
    }

    template<class DenseMatrixPolicy, class StatePolicy>
    static void CpuCalculateRateConstantsVectorized(const ReactionRateConstantStore& store, StatePolicy& state)
    {
      auto& v_rc = state.rate_constants_.AsVector();
      const auto& v_cp = state.custom_rate_parameters_.AsVector();
      constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
      const std::size_t n_cp = state.custom_rate_parameters_.NumColumns();
      const std::size_t n_analytic = store.lambda_offset();

      std::vector<double> cp(n_cp);
      std::vector<double> buf(n_analytic);

      for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
      {
        const std::size_t rc_base = i_group * state.rate_constants_.GroupSize();
        const std::size_t cp_base = i_group * state.custom_rate_parameters_.GroupSize();
        const std::size_t n_local = std::min(L, state.rate_constants_.NumRows() - i_group * L);

        for (std::size_t i_cell = 0; i_cell < n_local; ++i_cell)
        {
          const auto& cond = state.conditions_[i_group * L + i_cell];

          // Gather per-cell custom params from interleaved storage
          for (std::size_t j = 0; j < n_cp; ++j)
            cp[j] = v_cp[cp_base + j * L + i_cell];

          CpuCalculateOneCellRateConstants(store, cond.temperature_, cond.pressure_, cond.air_density_, cp.data(), buf.data());

          // Scatter results to interleaved layout
          for (std::size_t k = 0; k < n_analytic; ++k)
            v_rc[rc_base + k * L + i_cell] = buf[k];

          for (const auto& mult : store.parameterized_multipliers_)
            v_rc[rc_base + mult.rc_index * L + i_cell] *= mult.evaluate(cond);
        }
      }
    }
  };

}  // namespace micm
