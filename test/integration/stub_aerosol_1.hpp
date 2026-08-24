// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_1.hpp
/// @brief Stub aerosol model for integration testing
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Rate constants used by the test
constexpr micm::Real STUB1_RATE_CONSTANT_FO2_CORGE = 1e-3;
constexpr micm::Real STUB1_RATE_CONSTANT_BAZ_QUUX = 2e-3;

// First stubbed aerosol model implementation
//
// Mimics a single-moment two-mode aerosol model. Mode 1 contains only the first phase;
// mode 2 contains both phases. Two processes: FO2 gas → FO2 mode2 CORGE, and BAZ mode1 → BAZ mode2.
class StubAerosolModel
{
 public:
  struct RateConstants
  {
    micm::Real fo2_gas_to_mode2_corge_;
    micm::Real baz_mode1_to_mode2_quux_;
  };

  StubAerosolModel() = delete;
  StubAerosolModel(std::string name, const std::vector<micm::Phase>& phases, const RateConstants& rate_constants)
      : name_(std::move(name)),
        phases_(phases),
        rate_constants_(rate_constants)
  {
  }

  // State definition

  std::tuple<micm::Index, micm::Index> StateSize() const
  {
    EXPECT_EQ(phases_.size(), 2);
    micm::Index size = 0;
    size += phases_[0].StateSize();
    size += phases_[0].StateSize();
    size += phases_[1].StateSize();
    return { size, 0 };
  }

  std::set<std::string> StateVariableNames() const
  {
    std::set<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    for (const auto& name : phase1_names)
    {
      names.insert(name_ + ".MODE1." + name);
    }
    for (const auto& name : phase1_names)
    {
      names.insert(name_ + ".MODE2." + name);
    }
    for (const auto& name : phase2_names)
    {
      names.insert(name_ + ".MODE2." + name);
    }
    return names;
  }

  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  std::string Species(const micm::Index mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }

  // Process definition

  std::set<std::string> SpeciesUsed() const
  {
    return StateVariableNames();
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    std::set<std::pair<micm::Index, micm::Index>> elements;
    auto fo2_gas_it = state_indices.find("FO2");
    auto fo2_mode2_it = state_indices.find("STUB1.MODE2.CORGE.FO2");
    if (fo2_gas_it != state_indices.end() && fo2_mode2_it != state_indices.end())
    {
      elements.insert({ fo2_gas_it->second, fo2_gas_it->second });
      elements.insert({ fo2_mode2_it->second, fo2_gas_it->second });
    }
    auto baz_mode1_it = state_indices.find("STUB1.MODE1.QUUX.BAZ");
    auto baz_mode2_it = state_indices.find("STUB1.MODE2.QUUX.BAZ");
    if (baz_mode1_it != state_indices.end() && baz_mode2_it != state_indices.end())
    {
      elements.insert({ baz_mode1_it->second, baz_mode1_it->second });
      elements.insert({ baz_mode2_it->second, baz_mode1_it->second });
    }
    return elements;
  }

  /// Cache reactant/product state indices and pre-compute Jacobian flat indices for block 0.
  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    forcing_info_.clear();
    auto add = [&](const std::string& r, const std::string& p, micm::Real k)
    {
      auto ri = state_variable_indices.find(r);
      auto pi = state_variable_indices.find(p);
      if (ri != state_variable_indices.end() && pi != state_variable_indices.end())
      {
        forcing_info_.push_back({ ri->second,
                                  pi->second,
                                  k,
                                  jacobian.VectorIndex(0, ri->second, ri->second),
                                  jacobian.VectorIndex(0, pi->second, ri->second) });
      }
    };
    add("FO2", "STUB1.MODE2.CORGE.FO2", rate_constants_.fo2_gas_to_mode2_corge_);
    add("STUB1.MODE1.QUUX.BAZ", "STUB1.MODE2.QUUX.BAZ", rate_constants_.baz_mode1_to_mode2_quux_);
  }

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& /*conditions*/,
      DenseMatrixPolicy& /*state_parameters*/) const
  {
  }

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    for (const auto& info : forcing_info_)
    {
      const micm::Index reactant_idx = info.reactant_index_;
      const micm::Index product_idx = info.product_index_;
      const micm::Real k = info.rate_constant_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& forcing_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view) {
            forcing_view.ForEachRow(
                [k](micm::Real& f_reactant, micm::Real& f_product, const micm::Real& reactant)
                {
                  micm::Real rate = k * reactant;
                  f_reactant -= rate;
                  f_product += rate;
                },
                forcing_view.GetColumnView(reactant_idx),
                forcing_view.GetColumnView(product_idx),
                state_view.GetConstColumnView(reactant_idx));
          },
          forcing,
          state_variables)(forcing, state_variables);
    }
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractJacobianTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    for (const auto& info : forcing_info_)
    {
      const micm::Real k = info.rate_constant_;
      const micm::Index reactant_flat = info.reactant_reactant_flat_;
      const micm::Index product_flat = info.product_reactant_flat_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
            jacobian_view.ForEachBlock(
                [k](micm::Real& j_reactant, micm::Real& j_product)
                {
                  j_reactant -= -k;
                  j_product -= k;
                },
                jacobian_view.GetBlockView(reactant_flat),
                jacobian_view.GetBlockView(product_flat));
          },
          jacobian)(jacobian);
    }
    (void)state_variables;
  }

 private:
  struct ForcingInfo
  {
    micm::Index reactant_index_;
    micm::Index product_index_;
    micm::Real rate_constant_;
    micm::Index reactant_reactant_flat_;  // Jacobian flat id for (reactant, reactant)
    micm::Index product_reactant_flat_;   // Jacobian flat id for (product,  reactant)
  };

  std::string name_;
  std::vector<micm::Phase> phases_;
  RateConstants rate_constants_;
  std::vector<ForcingInfo> forcing_info_;
};
