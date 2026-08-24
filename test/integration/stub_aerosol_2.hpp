// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_2.hpp
/// @brief Stub aerosol model for integration testing
#pragma once

#include <micm/system/conditions.hpp>
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

// Second stubbed aerosol model implementation
//
// Mimics a two-moment three-mode aerosol model with temperature-dependent state parameters.
class AnotherStubAerosolModel
{
 public:
  AnotherStubAerosolModel() = delete;
  AnotherStubAerosolModel(std::string name, const std::vector<micm::Phase>& phases)
      : name_(std::move(name)),
        phases_(phases)
  {
  }

  // State definition

  std::tuple<micm::Index, micm::Index> StateSize() const
  {
    EXPECT_EQ(phases_.size(), 2);
    micm::Index size = 0;
    size += 1;
    size += phases_[0].StateSize();
    size += 1;
    size += phases_[1].StateSize();
    size += 1;
    size += phases_[0].StateSize();
    size += phases_[1].StateSize();
    return { size, 2 };
  }

  std::set<std::string> StateVariableNames() const
  {
    std::set<std::string> names;
    EXPECT_EQ(phases_.size(), 2);
    auto phase1_names = phases_[0].UniqueNames();
    auto phase2_names = phases_[1].UniqueNames();
    names.insert(name_ + ".MODE1.NUMBER");
    for (const auto& n : phase1_names)
    {
      names.insert(name_ + ".MODE1." + n);
    }
    names.insert(name_ + ".MODE2.NUMBER");
    for (const auto& n : phase2_names)
    {
      names.insert(name_ + ".MODE2." + n);
    }
    names.insert(name_ + ".MODE3.NUMBER");
    for (const auto& n : phase1_names)
    {
      names.insert(name_ + ".MODE3." + n);
    }
    for (const auto& n : phase2_names)
    {
      names.insert(name_ + ".MODE3." + n);
    }
    return names;
  }

  std::set<std::string> StateParameterNames() const
  {
    std::set<std::string> names;
    names.insert(name_ + ".PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
    names.insert(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
    return names;
  }

  std::string Species(const micm::Index mode, const micm::Phase& phase, const micm::Species& species) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + "." + phase.name_ + "." + species.name_;
  }

  std::string Number(const micm::Index mode) const
  {
    return name_ + ".MODE" + std::to_string(mode + 1) + ".NUMBER";
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
    auto fo2_mode2_it = state_indices.find("STUB2.MODE2.CORGE.FO2");
    auto baz_mode2_it = state_indices.find("STUB2.MODE2.CORGE.BAZ");
    if (fo2_mode2_it != state_indices.end() && baz_mode2_it != state_indices.end())
    {
      elements.insert({ fo2_mode2_it->second, fo2_mode2_it->second });
      elements.insert({ baz_mode2_it->second, fo2_mode2_it->second });
    }
    auto baz_mode3_it = state_indices.find("STUB2.MODE3.QUUX.BAZ");
    auto qux_mode3_it = state_indices.find("STUB2.MODE3.QUUX.QUX");
    if (baz_mode3_it != state_indices.end() && qux_mode3_it != state_indices.end())
    {
      elements.insert({ baz_mode3_it->second, baz_mode3_it->second });
      elements.insert({ qux_mode3_it->second, baz_mode3_it->second });
    }
    return elements;
  }

  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& state_parameter_indices,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    forcing_info_.clear();
    auto add = [&](const std::string& r, const std::string& p, const std::string& param)
    {
      auto ri = state_variable_indices.find(r);
      auto pi = state_variable_indices.find(p);
      auto pri = state_parameter_indices.find(param);
      if (ri != state_variable_indices.end() && pi != state_variable_indices.end() &&
          pri != state_parameter_indices.end())
      {
        forcing_info_.push_back({ ri->second,
                                  pi->second,
                                  pri->second,
                                  jacobian.VectorIndex(0, ri->second, ri->second),
                                  jacobian.VectorIndex(0, pi->second, ri->second) });
      }
    };
    add("STUB2.MODE2.CORGE.FO2", "STUB2.MODE2.CORGE.BAZ", name_ + ".PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT");
    add("STUB2.MODE3.QUUX.BAZ", "STUB2.MODE3.QUUX.QUX", name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");

    baz_to_qux_param_index_ = -1;
    auto it = state_parameter_indices.find(name_ + ".PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT");
    if (it != state_parameter_indices.end())
    {
      baz_to_qux_param_index_ = it->second;
    }
  }

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
      DenseMatrixPolicy& state_parameters) const
  {
    if (baz_to_qux_param_index_ < 0)
    {
      return;
    }
    const micm::Index idx = baz_to_qux_param_index_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
            const typename DenseMatrixPolicy::ViewType& state_param_view) {
          state_param_view.ForEachRow(
              [](const micm::Conditions& cond, micm::Real& rate) { rate = 0.005 * cond.temperature_; },
              conditions_view,
              state_param_view.GetColumnView(idx));
        },
        conditions,
        state_parameters)(conditions, state_parameters);
  }

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    for (const auto& info : forcing_info_)
    {
      const micm::Index reactant_idx = info.reactant_index_;
      const micm::Index product_idx = info.product_index_;
      const micm::Index param_idx = info.rate_param_index_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::ViewType& forcing_view,
              const typename DenseMatrixPolicy::ConstViewType& state_view,
              const typename DenseMatrixPolicy::ConstViewType& param_view) {
            forcing_view.ForEachRow(
                [](micm::Real& f_reactant,
                   micm::Real& f_product,
                   const micm::Real& reactant,
                   const micm::Real& rate_constant)
                {
                  micm::Real rate = rate_constant * reactant;
                  f_reactant -= rate;
                  f_product += rate;
                },
                forcing_view.GetColumnView(reactant_idx),
                forcing_view.GetColumnView(product_idx),
                state_view.GetConstColumnView(reactant_idx),
                param_view.GetConstColumnView(param_idx));
          },
          forcing,
          state_variables,
          state_parameters)(forcing, state_variables, state_parameters);
    }
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractJacobianTerms(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& /*state_variables*/,
      SparseMatrixPolicy& jacobian) const
  {
    for (const auto& info : forcing_info_)
    {
      const micm::Index param_idx = info.rate_param_index_;
      const micm::Index reactant_flat = info.reactant_reactant_flat_;
      const micm::Index product_flat = info.product_reactant_flat_;
      SparseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename SparseMatrixPolicy::ViewType& jacobian_view,
              const typename DenseMatrixPolicy::ConstViewType& param_view) {
            jacobian_view.ForEachBlock(
                [](micm::Real& j_reactant, micm::Real& j_product, const micm::Real& rate_constant)
                {
                  j_reactant -= -rate_constant;
                  j_product -= rate_constant;
                },
                jacobian_view.GetBlockView(reactant_flat),
                jacobian_view.GetBlockView(product_flat),
                param_view.GetConstColumnView(param_idx));
          },
          jacobian,
          state_parameters)(jacobian, state_parameters);
    }
  }

 private:
  struct ForcingInfo
  {
    micm::Index reactant_index_;
    micm::Index product_index_;
    micm::Index rate_param_index_;
    micm::Index reactant_reactant_flat_;
    micm::Index product_reactant_flat_;
  };

  std::string name_;
  std::vector<micm::Phase> phases_;
  std::vector<ForcingInfo> forcing_info_;
  int baz_to_qux_param_index_ = -1;
};
