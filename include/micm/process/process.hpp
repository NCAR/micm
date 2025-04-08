// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/error.hpp>

#include <memory>
#include <utility>
#include <vector>

enum class MicmProcessErrc
{
  TooManyReactantsForSurfaceReaction = MICM_PROCESS_ERROR_CODE_TOO_MANY_REACTANTS_FOR_SURFACE_REACTION
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmProcessErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmProcessErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_PROCESS;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmProcessErrc>(ev))
      {
        case MicmProcessErrc::TooManyReactantsForSurfaceReaction: return "A surface reaction can only have one reactant";
        default: return "Unknown error";
      }
    }
  };

  const MicmProcessErrorCategory MICM_PROCESS_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmProcessErrc e)
{
  return { static_cast<int>(e), MICM_PROCESS_ERROR };
}

namespace micm
{

  /// @brief An alias that allows making products easily
  using Yield = std::pair<micm::Species, double>;

  inline Yield Yields(const micm::Species& species, double yield)
  {
    return Yield(species, yield);
  };

  class ProcessBuilder;

  struct Process
  {
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;
    Phase phase_;

    /// @brief Recalculate the rate constants for each process for the current state
    /// @param processes The set of processes being solved
    /// @param state The solver state to update
    template<
        class DenseMatrixPolicy,
        class SparseMatrixPolicy,
        class LuDecompositionPolicy,
        class LMatrixPolicy,
        class UMatrixPolicy>
      requires(!VectorizableDense<DenseMatrixPolicy>)
    static void CalculateRateConstants(
        const std::vector<Process>& processes,
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state);
    template<
        class DenseMatrixPolicy,
        class SparseMatrixPolicy,
        class LuDecompositionPolicy,
        class LMatrixPolicy,
        class UMatrixPolicy>
      requires(VectorizableDense<DenseMatrixPolicy>)
    static void CalculateRateConstants(
        const std::vector<Process>& processes,
        State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state);

    friend class ProcessBuilder;
    static ProcessBuilder Create();
    Process(ProcessBuilder& builder);
    Process(const Process& other);

    Process(
        const std::vector<Species>& reactants,
        const std::vector<Yield>& products,
        std::unique_ptr<RateConstant> rate_constant,
        const Phase& phase)
        : reactants_(reactants),
          products_(products),
          rate_constant_(std::move(rate_constant)),
          phase_(phase)
    {
      if (dynamic_cast<SurfaceRateConstant*>(rate_constant_.get()))
      {
        if (reactants_.size() > 1)
        {
          throw std::system_error(make_error_code(MicmProcessErrc::TooManyReactantsForSurfaceReaction), "");
        }
      }
    }

    Process& operator=(const Process& other)
    {
      reactants_ = other.reactants_;
      products_ = other.products_;
      rate_constant_ = other.rate_constant_->Clone();
      phase_ = other.phase_;

      return *this;
    }
  };

  class ProcessBuilder
  {
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;
    Phase phase_;
    friend struct Process;

   public:
    operator Process() const
    {
      return Process(*this);
    }
    ProcessBuilder& SetReactants(const std::vector<Species>& reactants);
    ProcessBuilder& SetProducts(const std::vector<Yield>& products);
    ProcessBuilder& SetRateConstant(const RateConstant& rate_constant);
    ProcessBuilder& SetPhase(const Phase& phase);
  };

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
    requires(!VectorizableDense<DenseMatrixPolicy>)
  void Process::CalculateRateConstants(
      const std::vector<Process>& processes,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state)
  {
    MICM_PROFILE_FUNCTION();

    for (std::size_t i{}; i < state.custom_rate_parameters_.NumRows(); ++i)
    {
      const std::vector<double> custom_parameters = state.custom_rate_parameters_[i];
      std::vector<double>::const_iterator custom_parameters_iter = custom_parameters.begin();
      std::size_t i_rate_constant = 0;
      for (auto& process : processes)
      {
        double fixed_reactants = 1.0;
        for (auto& reactant : process.reactants_)
          if (reactant.IsParameterized())
            fixed_reactants *= reactant.parameterize_(state.conditions_[i]);
        state.rate_constants_[i][(i_rate_constant++)] =
            process.rate_constant_->Calculate(state.conditions_[i], custom_parameters_iter) * fixed_reactants;
        custom_parameters_iter += process.rate_constant_->SizeCustomParameters();
      }
    }
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
    requires(VectorizableDense<DenseMatrixPolicy>)
  void Process::CalculateRateConstants(
      const std::vector<Process>& processes,
      State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>& state)
  {
    MICM_PROFILE_FUNCTION();

    const auto& v_custom_parameters = state.custom_rate_parameters_.AsVector();
    auto& v_rate_constants = state.rate_constants_.AsVector();
    constexpr std::size_t L = DenseMatrixPolicy::GroupVectorSize();
    // loop over all rows
    for (std::size_t i_group = 0; i_group < state.rate_constants_.NumberOfGroups(); ++i_group)
    {
      std::size_t offset_rc = i_group * state.rate_constants_.GroupSize();
      std::size_t offset_params = i_group * state.custom_rate_parameters_.GroupSize();
      std::size_t rate_const_size = std::min(L, state.rate_constants_.NumRows() - (i_group * L));
      for (auto& process : processes)
      {
        std::vector<double> params(process.rate_constant_->SizeCustomParameters());
        for (std::size_t i_cell{}; i_cell < rate_const_size; ++i_cell)
        {
          for (std::size_t i_param = 0; i_param < params.size(); ++i_param)
          {
            params[i_param] = v_custom_parameters[offset_params + i_param * L + i_cell];
          }
          std::vector<double>::const_iterator custom_parameters_iter = params.begin();
          double fixed_reactants = 1.0;
          for (auto& reactant : process.reactants_)
            if (reactant.IsParameterized())
              fixed_reactants *= reactant.parameterize_(state.conditions_[i_group * L + i_cell]);
          v_rate_constants[offset_rc + i_cell] =
              process.rate_constant_->Calculate(state.conditions_[i_group * L + i_cell], custom_parameters_iter) *
              fixed_reactants;
        }
        offset_params += params.size() * L;
        offset_rc += L;
      }
    }
  }

  inline ProcessBuilder Process::Create()
  {
    return ProcessBuilder{};
  };

  inline Process::Process(ProcessBuilder& builder)
      : Process(builder.reactants_, builder.products_, std::move(builder.rate_constant_), builder.phase_)
  {
  }

  inline Process::Process(const Process& other)
      : reactants_(other.reactants_),
        products_(other.products_),
        rate_constant_(other.rate_constant_ ? other.rate_constant_->Clone() : nullptr),
        phase_(other.phase_)
  {
  }

  inline ProcessBuilder& ProcessBuilder::SetReactants(const std::vector<Species>& reactants)
  {
    reactants_ = reactants;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::SetProducts(const std::vector<Yield>& products)
  {
    products_ = products;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::SetRateConstant(const RateConstant& rate_constant)
  {
    rate_constant_ = rate_constant.Clone();
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::SetPhase(const Phase& phase)
  {
    phase_ = phase;
    return *this;
  }
}  // namespace micm
