// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/external_model.hpp>
#include <micm/util/types.hpp>

#include <array>
#include <memory>
#include <tuple>
#include <utility>

namespace micm
{
  /// @brief Wraps an inner rates policy and a shared tuple of concrete external models.
  ///
  /// Solve-time methods first delegate to the inner (built-in) rates policy, then dispatch
  /// directly on each external model that satisfies `HasProcesses`.
  template<class InnerRates, class... ExternalModels>
  class RatesBundle
  {
   public:
    using ModelsTuple = std::tuple<ExternalModels...>;

    RatesBundle() = default;

    RatesBundle(InnerRates&& inner, std::shared_ptr<ModelsTuple> models)
        : inner_(std::move(inner)),
          models_(std::move(models))
    {
    }

    RatesBundle(RatesBundle&&) noexcept = default;
    RatesBundle& operator=(RatesBundle&&) noexcept = default;
    RatesBundle(const RatesBundle&) = delete;
    RatesBundle& operator=(const RatesBundle&) = delete;

    InnerRates& inner()
    {
      return inner_;
    }
    const InnerRates& inner() const
    {
      return inner_;
    }

    template<class State, class DenseMatrixPolicy>
    void AddForcingTerms(const State& state, const DenseMatrixPolicy& Y, DenseMatrixPolicy& forcing) const
    {
      inner_.AddForcingTerms(state, Y, forcing);
      InvokeProcesses([&](const auto& m)
                      { m.AddForcingTerms(state.custom_rate_parameters_, Y, forcing); });
    }

    template<class State, class DenseMatrixPolicy, class SparseMatrixPolicy>
    void SubtractJacobianTerms(const State& state, const DenseMatrixPolicy& Y, SparseMatrixPolicy& jacobian) const
    {
      inner_.SubtractJacobianTerms(state, Y, jacobian);
      InvokeProcesses([&](const auto& m)
                      { m.SubtractJacobianTerms(state.custom_rate_parameters_, Y, jacobian); });
    }

    /// @brief Called before each solve to refresh temperature-/pressure-dependent parameters.
    template<class ConditionsVector, class DenseMatrixPolicy>
    void UpdateStateParameters(const ConditionsVector& conditions, DenseMatrixPolicy& state_parameters) const
    {
      InvokeProcesses([&](const auto& m) { m.UpdateStateParameters(conditions, state_parameters); });
    }

    /// @brief Forward CUDA store upload to the inner rates policy when it supports it.
    template<class Store>
    void BuildCudaStore(const Store& store)
      requires requires(InnerRates& r) { r.BuildCudaStore(store); }
    {
      inner_.BuildCudaStore(store);
    }

    template<class Store, class State>
    void GpuCalculateRateConstants(const Store& store, State& state)
      requires requires(InnerRates& r) { r.GpuCalculateRateConstants(store, state); }
    {
      inner_.GpuCalculateRateConstants(store, state);
    }

   private:
    template<class F>
    void InvokeProcesses(F&& f) const
    {
      if (!models_)
      {
        return;
      }
      std::apply(
          [&](const auto&... m)
          {
            auto invoke = [&](const auto& model)
            {
              using M = std::decay_t<decltype(model)>;
              if constexpr (HasProcesses<M>)
              {
                f(model);
              }
            };
            (invoke(m), ...);
          },
          *models_);
    }

    InnerRates inner_{};
    std::shared_ptr<ModelsTuple> models_;
  };

  /// @brief Wraps an inner constraint set and a shared tuple of concrete external models.
  ///
  /// Only models that both satisfy `HasConstraints` AND report a non-empty algebraic-variable
  /// set at build time contribute at solve time. The `active_` mask is populated by the builder
  /// so runtime-configurable models (that opt out via empty names) are cheaply skipped.
  template<class InnerConstraints, class... ExternalModels>
  class ConstraintBundle
  {
   public:
    static constexpr std::size_t model_count = sizeof...(ExternalModels);
    using ModelsTuple = std::tuple<ExternalModels...>;
    using ActiveMask = std::array<bool, model_count>;

    ConstraintBundle() = default;

    ConstraintBundle(InnerConstraints&& inner, std::shared_ptr<ModelsTuple> models, ActiveMask active)
        : inner_(std::move(inner)),
          models_(std::move(models)),
          active_(active)
    {
    }

    ConstraintBundle(ConstraintBundle&&) noexcept = default;
    ConstraintBundle& operator=(ConstraintBundle&&) noexcept = default;
    ConstraintBundle(const ConstraintBundle&) = default;
    ConstraintBundle& operator=(const ConstraintBundle&) = default;

    InnerConstraints& inner()
    {
      return inner_;
    }
    const InnerConstraints& inner() const
    {
      return inner_;
    }

    /// Number of algebraic constraint rows (built-in + active external models).
    Index Size() const
    {
      return static_cast<Index>(inner_.AlgebraicVariableIds().size());
    }

    template<class DenseMatrixPolicy>
    void AddForcingTerms(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        DenseMatrixPolicy& forcing) const
    {
      inner_.AddForcingTerms(state_variables, state_parameters, forcing);
      InvokeConstraints([&](const auto& m) { m.AddConstraintResidual(state_parameters, state_variables, forcing); });
    }

    template<class DenseMatrixPolicy, class SparseMatrixPolicy>
    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_variables,
        const DenseMatrixPolicy& state_parameters,
        SparseMatrixPolicy& jacobian) const
    {
      inner_.SubtractJacobianTerms(state_variables, state_parameters, jacobian);
      InvokeConstraints([&](const auto& m)
                        { m.SubtractConstraintJacobian(state_parameters, state_variables, jacobian); });
    }

    template<class DenseMatrixPolicy>
    void SetAlgebraicErrors(DenseMatrixPolicy& Yerror, const DenseMatrixPolicy& Y, const DenseMatrixPolicy& Ynew) const
    {
      inner_.SetAlgebraicErrors(Yerror, Y, Ynew);
    }

    template<class ConditionsVector, class DenseMatrixPolicy>
    void UpdateStateParameters(const ConditionsVector& conditions, DenseMatrixPolicy& state_parameters) const
    {
      inner_.UpdateStateParameters(conditions, state_parameters);
      InvokeConstraints([&](const auto& m) { m.UpdateConstraintStateParameters(conditions, state_parameters); });
    }

    /// @brief Diagnose constraint parameters from state at the start of each Solve().
    template<class DenseMatrixPolicy>
    void InitializeConstraintParameters(const DenseMatrixPolicy& state_variables, DenseMatrixPolicy& state_parameters) const
    {
      InvokeConstraints(
          [&](const auto& m)
          {
            using M = std::decay_t<decltype(m)>;
            if constexpr (HasInitializeConstraintParameters<M>)
            {
              m.InitializeConstraintParameters(state_variables, state_parameters);
            }
          });
    }

   private:
    template<class F>
    void InvokeConstraints(F&& f) const
    {
      if (!models_)
      {
        return;
      }
      InvokeConstraintsImpl(std::forward<F>(f), std::make_index_sequence<model_count>{});
    }

    template<class F, std::size_t... I>
    void InvokeConstraintsImpl(F&& f, std::index_sequence<I...>) const
    {
      (InvokeOne<I>(f), ...);
    }

    template<std::size_t I, class F>
    void InvokeOne(F& f) const
    {
      using M = std::tuple_element_t<I, ModelsTuple>;
      if constexpr (HasConstraints<M>)
      {
        if (active_[I])
        {
          f(std::get<I>(*models_));
        }
      }
    }

    InnerConstraints inner_{};
    std::shared_ptr<ModelsTuple> models_;
    ActiveMask active_{};
  };
}  // namespace micm
