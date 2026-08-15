// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file reaction_rate_store.hpp
/// @brief Structure-of-arrays store for all reaction rate constant parameters.
///
/// Reactions must be stable-sorted by RateConstantTypeOrder before BuildFrom is called
/// (the SolverBuilder does this).  Each type occupies a contiguous block; offset helpers
/// below give the start of each block within state.rate_constants_[cell].

#define _USE_MATH_DEFINES

#include <micm/process/process.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/solver/state.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/micm_exception.hpp>
#include <micm/util/parameterized_function.hpp>
#include <micm/util/property_keys.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
  /// Templated on DenseMatrixPolicy to pick up VectorType<T> — for Kokkos builds
  /// this is KokkosPaddedVector<T> which has CopyToDevice()/GetView() device support,
  /// matching the same Vector/VectorView pattern used in the LU decomposers.
  template<class DenseMatrixPolicy>
  struct ReactionRateConstantStore
  {
    template<class T>
    using Vector = typename DenseMatrixPolicy::template VectorType<T>;
    template<class T>
    using VectorView = typename Vector<T>::ConstViewType;

    // ----------------------------------------------------------------
    // Parameterized-species multipliers (device-safe POD)
    // ----------------------------------------------------------------

    /// @brief One entry per reaction with at least one parameterized reactant.
    ///        Trivially copyable so it can live in a device Kokkos::View.
    struct ParameterizedMultiplier
    {
      /// @brief Maximum number of parameterized reactants supported in a single
      ///        reaction.  If a reaction exceeds this, BuildFrom throws.
      static constexpr Index kMaxFuncs = 4;

      ParameterizedFunction funcs_[kMaxFuncs]{};
      Index n_funcs_ = 0;
      Index rc_index_ = 0;

      MICM_DEVICE_FUNCTION Real Evaluate(const Conditions& c) const
      {
        Real v = 1.0;
        for (Index i = 0; i < n_funcs_; ++i)
        {
          v *= funcs_[i](c);
        }
        return v;
      }
    };

    // ----------------------------------------------------------------
    // Analytic types (GPU-safe parameter structs) — stored in device-capable vectors
    // ----------------------------------------------------------------
    Vector<ArrheniusRateConstantParameters> arrhenius_;
    Vector<TroeRateConstantParameters> troe_;
    Vector<TernaryChemicalActivationRateConstantParameters> ternary_;
    Vector<BranchedRateConstantParameters> branched_;
    Vector<TunnelingRateConstantParameters> tunneling_;
    Vector<TaylorSeriesRateConstantParameters> taylor_;
    Vector<ReversibleRateConstantParameters> reversible_;
    Vector<UserDefinedRateConstantData> user_defined_;
    Vector<SurfaceRateConstantData> surface_;
    Vector<ParameterizedMultiplier> parameterized_multipliers_;

    // ----------------------------------------------------------------
    // Device-accessible views (rebuilt in CopyToDevice after population)
    // ----------------------------------------------------------------
    struct Views
    {
      VectorView<ArrheniusRateConstantParameters> arrhenius_;
      VectorView<TroeRateConstantParameters> troe_;
      VectorView<TernaryChemicalActivationRateConstantParameters> ternary_;
      VectorView<BranchedRateConstantParameters> branched_;
      VectorView<TunnelingRateConstantParameters> tunneling_;
      VectorView<TaylorSeriesRateConstantParameters> taylor_;
      VectorView<ReversibleRateConstantParameters> reversible_;
      VectorView<UserDefinedRateConstantData> user_defined_;
      VectorView<SurfaceRateConstantData> surface_;
      VectorView<ParameterizedMultiplier> parameterized_multipliers_;

      Views() = default;

      Views(
          const Vector<ArrheniusRateConstantParameters>& arr,
          const Vector<TroeRateConstantParameters>& troe,
          const Vector<TernaryChemicalActivationRateConstantParameters>& tern,
          const Vector<BranchedRateConstantParameters>& bran,
          const Vector<TunnelingRateConstantParameters>& tunn,
          const Vector<TaylorSeriesRateConstantParameters>& tayl,
          const Vector<ReversibleRateConstantParameters>& rev,
          const Vector<UserDefinedRateConstantData>& ud,
          const Vector<SurfaceRateConstantData>& surf,
          const Vector<ParameterizedMultiplier>& pmult)
          : arrhenius_(arr.GetView()),
            troe_(troe.GetView()),
            ternary_(tern.GetView()),
            branched_(bran.GetView()),
            tunneling_(tunn.GetView()),
            taylor_(tayl.GetView()),
            reversible_(rev.GetView()),
            user_defined_(ud.GetView()),
            surface_(surf.GetView()),
            parameterized_multipliers_(pmult.GetView())
      {
      }
    };
    Views views_;

    void CopyToDevice()
    {
      arrhenius_.CopyToDevice();
      troe_.CopyToDevice();
      ternary_.CopyToDevice();
      branched_.CopyToDevice();
      tunneling_.CopyToDevice();
      taylor_.CopyToDevice();
      reversible_.CopyToDevice();
      user_defined_.CopyToDevice();
      surface_.CopyToDevice();
      parameterized_multipliers_.CopyToDevice();
      views_ = Views(
          arrhenius_,
          troe_,
          ternary_,
          branched_,
          tunneling_,
          taylor_,
          reversible_,
          user_defined_,
          surface_,
          parameterized_multipliers_);
    }

    // ----------------------------------------------------------------
    // Precomputed contiguous-block offsets into state.rate_constants_[cell]
    // ----------------------------------------------------------------
    Index off_troe_{ 0 }, off_tern_{ 0 }, off_bran_{ 0 }, off_tunn_{ 0 }, off_tayl_{ 0 }, off_rev_{ 0 }, off_ud_{ 0 },
        off_surf_{ 0 }, off_lambda_{ 0 };

    Index TroeOffset() const
    {
      return off_troe_;
    }
    Index TernaryOffset() const
    {
      return off_tern_;
    }
    Index BranchedOffset() const
    {
      return off_bran_;
    }
    Index TunnelingOffset() const
    {
      return off_tunn_;
    }
    Index TaylorOffset() const
    {
      return off_tayl_;
    }
    Index ReversibleOffset() const
    {
      return off_rev_;
    }
    Index UserDefinedOffset() const
    {
      return off_ud_;
    }
    Index SurfaceOffset() const
    {
      return off_surf_;
    }
    Index LambdaOffset() const
    {
      return off_lambda_;
    }

    // ----------------------------------------------------------------
    // CPU-only lambda entries
    // ----------------------------------------------------------------

    struct LambdaEntry
    {
      LambdaRateConstantParameters* source_;  ///< Non-owning; valid for the lifetime of the owning Solver.
      Index rc_index_;                        ///< Column index in state.rate_constants_[cell].
    };
    std::vector<LambdaEntry> lambda_entries_;

    // ================================================================
    // Factory
    // ================================================================

    /// @brief Build a ReactionRateConstantStore from a sorted process list.
    /// @param processes  Non-const ref so LambdaRateConstantParameters pointers remain mutable at runtime.
    static ReactionRateConstantStore<DenseMatrixPolicy> BuildFrom(std::vector<Process>& processes)
    {
      ReactionRateConstantStore<DenseMatrixPolicy> store;
      // Accumulate into plain std::vector (push_back friendly), then assign to Vector<T>
      std::vector<ArrheniusRateConstantParameters> arrhenius_tmp;
      std::vector<TroeRateConstantParameters> troe_tmp;
      std::vector<TernaryChemicalActivationRateConstantParameters> ternary_tmp;
      std::vector<BranchedRateConstantParameters> branched_tmp;
      std::vector<TunnelingRateConstantParameters> tunneling_tmp;
      std::vector<TaylorSeriesRateConstantParameters> taylor_tmp;
      std::vector<ReversibleRateConstantParameters> reversible_tmp;
      std::vector<UserDefinedRateConstantData> user_defined_tmp;
      std::vector<SurfaceRateConstantData> surface_tmp;
      std::vector<ParameterizedMultiplier> parameterized_multipliers_tmp;

      Index rc_index = 0;
      Index custom_param_off = 0;

      for (auto& process : processes)
      {
        auto& reaction = process.process_;

        {  // parameterized-reactant multiplier
          ParameterizedMultiplier mult{};
          mult.rc_index_ = rc_index;
          for (const auto& reactant : reaction.reactants_)
          {
            if (reactant.IsParameterized())
            {
              if (mult.n_funcs_ >= ParameterizedMultiplier::kMaxFuncs)
              {
                throw MicmException(
                    MICM_ERROR_CATEGORY_SPECIES,
                    MICM_SPECIES_ERROR_CODE_PROPERTY_NOT_FOUND,
                    "Reaction has more parameterized reactants than ParameterizedMultiplier::kMaxFuncs; "
                    "increase kMaxFuncs.");
              }
              mult.funcs_[mult.n_funcs_++] = reactant.parameterize_;
            }
          }

          if (mult.n_funcs_ > 0)
          {
            parameterized_multipliers_tmp.push_back(mult);
          }
        }

        Index n_custom = 0;

        RateConstantVariant& rc = reaction.rate_constant_;

        if (auto* p = std::get_if<ArrheniusRateConstantParameters>(&rc))
        {
          arrhenius_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<TroeRateConstantParameters>(&rc))
        {
          troe_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<TernaryChemicalActivationRateConstantParameters>(&rc))
        {
          ternary_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<BranchedRateConstantParameters>(&rc))
        {
          BranchedRateConstantParameters params = *p;
          // Pre-compute derived fields needed by CalculateBranched
          params.k0_ = 2.0e-22 * constants::AVOGADRO_CONSTANT * 1.0e-6 * std::exp(static_cast<Real>(params.n_));
          Real air_ref = 2.45e19 / constants::AVOGADRO_CONSTANT * 1.0e6;
          Real a = params.k0_ * air_ref;
          Real b = 0.43 * std::pow(293.0 / 298.0, -8.0);
          Real A_val = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
          params.z_ = A_val * (1.0 - params.a0_) / params.a0_;
          branched_tmp.push_back(params);
        }
        else if (auto* p = std::get_if<TunnelingRateConstantParameters>(&rc))
        {
          tunneling_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<TaylorSeriesRateConstantParameters>(&rc))
        {
          taylor_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<ReversibleRateConstantParameters>(&rc))
        {
          reversible_tmp.push_back(*p);
        }
        else if (auto* p = std::get_if<UserDefinedRateConstantParameters>(&rc))
        {
          UserDefinedRateConstantData data;
          data.scaling_factor_ = p->scaling_factor_;
          data.custom_param_index_ = custom_param_off;
          user_defined_tmp.push_back(data);
          n_custom = 1;
        }
        else if (auto* p = std::get_if<SurfaceRateConstantParameters>(&rc))
        {
          SurfaceRateConstantData data;
          if (!p->phase_species_.diffusion_coefficient_.has_value())
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_SPECIES,
                MICM_SPECIES_ERROR_CODE_PROPERTY_NOT_FOUND,
                "Diffusion coefficient for species '" + p->phase_species_.species_.name_ + "' is not defined");
          }
          data.diffusion_coefficient_ = p->phase_species_.diffusion_coefficient_.value();
          auto mw = p->phase_species_.species_.GetProperty<Real>(property_keys::MOLECULAR_WEIGHT);
          data.mean_free_speed_factor_ = 8.0 * constants::GAS_CONSTANT / (M_PI * mw);
          data.reaction_probability_ = p->reaction_probability_;
          data.custom_param_base_index_ = custom_param_off;
          surface_tmp.push_back(data);
          n_custom = 2;
        }
        else if (auto* p = std::get_if<LambdaRateConstantParameters>(&rc))
        {
          store.lambda_entries_.push_back({ p, rc_index });
        }

        custom_param_off += n_custom;
        ++rc_index;
      }

      // Assign temp vectors to device-capable Vector<T> members and upload
      store.arrhenius_ = Vector<ArrheniusRateConstantParameters>(arrhenius_tmp);
      store.troe_ = Vector<TroeRateConstantParameters>(troe_tmp);
      store.ternary_ = Vector<TernaryChemicalActivationRateConstantParameters>(ternary_tmp);
      store.branched_ = Vector<BranchedRateConstantParameters>(branched_tmp);
      store.tunneling_ = Vector<TunnelingRateConstantParameters>(tunneling_tmp);
      store.taylor_ = Vector<TaylorSeriesRateConstantParameters>(taylor_tmp);
      store.reversible_ = Vector<ReversibleRateConstantParameters>(reversible_tmp);
      store.user_defined_ = Vector<UserDefinedRateConstantData>(user_defined_tmp);
      store.surface_ = Vector<SurfaceRateConstantData>(surface_tmp);
      store.parameterized_multipliers_ = Vector<ParameterizedMultiplier>(parameterized_multipliers_tmp);

      // Precompute offsets
      store.off_troe_ = static_cast<Index>(arrhenius_tmp.size());
      store.off_tern_ = store.off_troe_ + static_cast<Index>(troe_tmp.size());
      store.off_bran_ = store.off_tern_ + static_cast<Index>(ternary_tmp.size());
      store.off_tunn_ = store.off_bran_ + static_cast<Index>(branched_tmp.size());
      store.off_tayl_ = store.off_tunn_ + static_cast<Index>(tunneling_tmp.size());
      store.off_rev_ = store.off_tayl_ + static_cast<Index>(taylor_tmp.size());
      store.off_ud_ = store.off_rev_ + static_cast<Index>(reversible_tmp.size());
      store.off_surf_ = store.off_ud_ + static_cast<Index>(user_defined_tmp.size());
      store.off_lambda_ = store.off_surf_ + static_cast<Index>(surface_tmp.size());

      store.CopyToDevice();
      return store;
    }

    /// @brief Evaluate all lambda rate constants into state.rate_constants_.
    ///        Called prior to calculating device-compatible rate constants
    template<class StatePolicy>
    static void CalculateCpuRateConstants(const ReactionRateConstantStore<DenseMatrixPolicy>& store, StatePolicy& state)
    {
      if (store.lambda_entries_.empty())
      {
        return;
      }

      DenseMatrixPolicy::HostFunction(
          [&store](auto&& rc_view, const auto& conditions)
          {
            for (const auto& entry : store.lambda_entries_)
            {
              rc_view.ForEachRow(
                  [&entry](Real& rate_constant, const Conditions& conditions)
                  { rate_constant = entry.source_->lambda_function_(conditions); },
                  rc_view.GetColumnView(entry.rc_index_),
                  conditions);
            }
          },
          state.rate_constants_,
          state.conditions_)(state.rate_constants_, state.conditions_);
      state.rate_constants_.CopyToDevice();
    }

    /// @brief Calculate all analytic rate constants into state.rate_constants_.
    ///        Lambda entries are untouched; parameterized multipliers applied last.
    ///        Uses DenseMatrixPolicy::Function so a single implementation handles both
    ///        Matrix (scalar) and VectorMatrix (interleaved) layouts. Each reaction type
    ///        is computed across all cells at once via ForEachRow, which is more
    ///        SIMD-friendly than the previous per-cell loop.
    template<class StatePolicy>
    static void CalculateRateConstants(const ReactionRateConstantStore<DenseMatrixPolicy>& store, StatePolicy& state)
    {
      using DenseMatrix = DenseMatrixPolicy;
      CalculateCpuRateConstants(store, state);
      // Capture device-accessible views and precomputed offsets/sizes for the GPU kernel.
      // views_ members are VectorView<T> (KokkosPaddedVector::ConstDeviceView for Kokkos builds)
      // which hold a Kokkos::View pointing to device-resident data uploaded by CopyToDevice().
      const auto& v = store.views_;
      const auto n_arr = static_cast<Index>(v.arrhenius_.size());
      const auto n_troe = static_cast<Index>(v.troe_.size());
      const auto n_tern = static_cast<Index>(v.ternary_.size());
      const auto n_bran = static_cast<Index>(v.branched_.size());
      const auto n_tunn = static_cast<Index>(v.tunneling_.size());
      const auto n_tayl = static_cast<Index>(v.taylor_.size());
      const auto n_rev = static_cast<Index>(v.reversible_.size());
      const auto n_ud = static_cast<Index>(v.user_defined_.size());
      const auto n_surf = static_cast<Index>(v.surface_.size());
      const auto n_mult = static_cast<Index>(v.parameterized_multipliers_.size());
      const Index off_troe = store.off_troe_, off_tern = store.off_tern_, off_bran = store.off_bran_,
                  off_tunn = store.off_tunn_, off_tayl = store.off_tayl_, off_rev = store.off_rev_, off_ud = store.off_ud_,
                  off_surf = store.off_surf_;
      DenseMatrix::Function(
          MICM_LAMBDA(
              typename DenseMatrix::ViewType rate_constants,
              typename DenseMatrix::ConstViewType parameters,
              typename Vector<Conditions>::ConstViewType conditions) {
            for (Index i = 0; i < n_arr; ++i)
            {
              const auto& p = v.arrhenius_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions)
                  { out = CalculateArrhenius(p, conditions.temperature_, conditions.pressure_); },
                  rate_constants.GetColumnView(i),
                  conditions);
            }
            for (Index i = 0; i < n_troe; ++i)
            {
              const auto& p = v.troe_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions)
                  { out = CalculateTroe(p, conditions.temperature_, conditions.air_density_); },
                  rate_constants.GetColumnView(off_troe + i),
                  conditions);
            }
            for (Index i = 0; i < n_tern; ++i)
            {
              const auto& p = v.ternary_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions)
                  { out = CalculateTernaryChemicalActivation(p, conditions.temperature_, conditions.air_density_); },
                  rate_constants.GetColumnView(off_tern + i),
                  conditions);
            }
            for (Index i = 0; i < n_bran; ++i)
            {
              const auto& p = v.branched_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions)
                  { out = CalculateBranched(p, conditions.temperature_, conditions.air_density_); },
                  rate_constants.GetColumnView(off_bran + i),
                  conditions);
            }
            for (Index i = 0; i < n_tunn; ++i)
            {
              const auto& p = v.tunneling_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions) { out = CalculateTunneling(p, conditions.temperature_); },
                  rate_constants.GetColumnView(off_tunn + i),
                  conditions);
            }
            for (Index i = 0; i < n_tayl; ++i)
            {
              const auto& p = v.taylor_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions)
                  { out = CalculateTaylorSeries(p, conditions.temperature_, conditions.pressure_); },
                  rate_constants.GetColumnView(off_tayl + i),
                  conditions);
            }
            for (Index i = 0; i < n_rev; ++i)
            {
              const auto& p = v.reversible_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions) { out = CalculateReversible(p, conditions.temperature_); },
                  rate_constants.GetColumnView(off_rev + i),
                  conditions);
            }
            for (Index i = 0; i < n_ud; ++i)
            {
              const auto& p = v.user_defined_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Real& parameter) { out = CalculateUserDefined(p, parameter); },
                  rate_constants.GetColumnView(off_ud + i),
                  parameters.GetConstColumnView(p.custom_param_index_));
            }
            for (Index i = 0; i < n_surf; ++i)
            {
              const auto& p = v.surface_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& conditions, const Real& radius, const Real& num_conc)
                  { out = CalculateSurfaceOne(p, conditions.temperature_, radius, num_conc); },
                  rate_constants.GetColumnView(off_surf + i),
                  conditions,
                  parameters.GetConstColumnView(p.custom_param_base_index_),
                  parameters.GetConstColumnView(p.custom_param_base_index_ + 1));
            }
            for (Index i = 0; i < n_mult; ++i)
            {
              const auto& mult = v.parameterized_multipliers_[i];
              rate_constants.ForEachRow(
                  [=](Real& out, const Conditions& c) { out *= mult.Evaluate(c); },
                  rate_constants.GetColumnView(mult.rc_index_),
                  conditions);
            }
          },
          state.rate_constants_,
          state.custom_rate_parameters_,
          state.conditions_)(state.rate_constants_, state.custom_rate_parameters_, state.conditions_);
    }
  };

}  // namespace micm
