/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <memory>
#include <vector>
#include <utility>
#include <micm/system/species.hpp>
#include <micm/system/phase.hpp>
#include <micm/process/rate_constant.hpp>

namespace micm
{

  using Yield = std::pair<micm::Species, double>;

  Yield yields(micm::Species species, double yield) {
    return Yield(species, yield);
  };

  class ProcessBuilder;

  struct Process {
    const std::vector<Species> reactants_;
    const std::vector<Yield> products_;
    const std::unique_ptr<RateConstant> rate_constant_;
    const Phase phase_;
    friend class ProcessBuilder;
    static ProcessBuilder create();
    Process(ProcessBuilder& builder);
    Process(const Process& other);
  };

  class ProcessBuilder {
    std::vector<Species> reactants_;
    std::vector<Yield> products_;
    std::unique_ptr<RateConstant> rate_constant_;
    Phase phase_;  
    friend struct Process;
  public:
    operator Process() const { return Process(*this); }
    ProcessBuilder& reactants(const std::vector<Species>& reactants);
    ProcessBuilder& products(const std::vector<Yield>& products);
    ProcessBuilder& rate_constant(const RateConstant& rate_constant);
    ProcessBuilder& phase(const Phase& phase);
  };

  inline ProcessBuilder Process::create() { return ProcessBuilder{}; };

  Process::Process(ProcessBuilder& builder)
    : reactants_(builder.reactants_)
    , products_(builder.products_)
    , rate_constant_(std::move(builder.rate_constant_))
    , phase_(builder.phase_)
    {}

  Process::Process(const Process& other)
    : reactants_(other.reactants_)
    , products_(other.products_)
    , rate_constant_(other.rate_constant_->clone())
    , phase_(other.phase_)
    {}

  inline ProcessBuilder& ProcessBuilder::reactants(const std::vector<Species>& reactants) {
    reactants_ = reactants;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::products(const std::vector<Yield>& products) {
    products_ = products;
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::rate_constant(const RateConstant& rate_constant) {
    rate_constant_ = rate_constant.clone();
    return *this;
  }

  inline ProcessBuilder& ProcessBuilder::phase(const Phase& phase) {
    phase_ = phase;
    return *this;
  }
}  // namespace micm
