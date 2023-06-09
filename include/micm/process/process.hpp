/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <memory>
#include <micm/process/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <utility>
#include <vector>

namespace micm
{

  using Yield = std::pair<micm::Species, double>;

  Yield yields(micm::Species species, double yield)
  {
    return Yield(species, yield);
  };

  class ProcessBuilder;
  
  struct Process
  {
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::shared_ptr<RateConstant> rate_constant_;
    Phase phase_;

    /// @brief Update the solver state rate constants
    /// @param processes The set of processes being solved
    /// @param state The solver state to update
    static void UpdateState(const std::vector<Process>& processes, State& state);

    friend class ProcessBuilder;
    static ProcessBuilder create();
    Process(ProcessBuilder& builder);
    Process(const Process& other);
  

    Process(const std::vector<Species>& reactants, const std::vector<Yield>& products, std::shared_ptr<RateConstant> rate_constant, const Phase& phase)
        : reactants_(reactants),
        products_(products),
        rate_constant_(rate_constant),  // TODO: jiwon 6/7 - need to be ptr, shared ptr okay?
        phase_(phase)
    {
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
    ProcessBuilder& reactants(const std::vector<Species>& reactants);
    ProcessBuilder& products(const std::vector<Yield>& products);
    ProcessBuilder& rate_constant(const RateConstant& rate_constant);
    ProcessBuilder& phase(const Phase& phase);
  };

  void Process::UpdateState(const std::vector<Process>& processes, State& state)
  {
    for (std::size_t i{}; i < state.custom_rate_parameters_.size(); ++i)
    {
      std::vector<double>::const_iterator custom_parameters = state.custom_rate_parameters_[i].begin();
      std::vector<double>::iterator rate_constant = state.rate_constants_[i].begin();
      for (auto& process : processes)
      {
        *(rate_constant++) = process.rate_constant_->calculate(state.conditions_[i], custom_parameters);
        custom_parameters += process.rate_constant_->SizeCustomParameters();
      }
    }
  }

  inline ProcessBuilder Process::create()
  {
    return ProcessBuilder{};
  };

  Process::Process(ProcessBuilder& builder)
      : reactants_(builder.reactants_),
        products_(builder.products_),
        rate_constant_(std::move(builder.rate_constant_)),
        phase_(builder.phase_)
  {
  }

  Process::Process(const Process& other)
      : reactants_(other.reactants_),
        products_(other.products_),
        rate_constant_(other.rate_constant_ ? other.rate_constant_->clone() : nullptr),
        phase_(other.phase_)
  {
  }

  inline ProcessBuilder& ProcessBuilder::reactants(const std::vector<Species>& reactants)
  {
    reactants_ = reactants;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::products(const std::vector<Yield>& products)
  {
    products_ = products;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::rate_constant(const RateConstant& rate_constant)
  {
    rate_constant_ = rate_constant.clone();
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::phase(const Phase& phase)
  {
    phase_ = phase;
    return *this;
  }
}  // namespace micm
