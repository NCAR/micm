// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_with_constraints.hpp
/// @brief Stub aerosol model with constraints for integration testing
#pragma once

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/types.hpp>

#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

/// A stub external model providing both processes and constraints.
///
/// A_GAS → A_AQ at rate k. When `total_mass` is provided, adds the mass-conservation
/// constraint [A_GAS] + [A_AQ] = total (DAE mode).
class StubAerosolWithConstraints
{
 public:
  StubAerosolWithConstraints() = delete;

  StubAerosolWithConstraints(micm::Real rate_constant, std::optional<micm::Real> total_mass = std::nullopt)
      : rate_constant_(rate_constant),
        total_mass_(total_mass)
  {
  }

  // State definition

  std::tuple<micm::Index, micm::Index> StateSize() const
  {
    return { 1, 0 };
  }
  std::set<std::string> StateVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }
  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  // Process definition

  std::set<std::string> SpeciesUsed() const
  {
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    std::set<std::pair<micm::Index, micm::Index>> elements;
    auto i_gas = state_indices.find("A_GAS");
    auto i_aq = state_indices.find("AEROSOL.A_AQ");
    if (i_gas != state_indices.end() && i_aq != state_indices.end())
    {
      elements.insert({ i_gas->second, i_gas->second });
      elements.insert({ i_aq->second, i_gas->second });
    }
    return elements;
  }

  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    i_gas_ = state_variable_indices.at("A_GAS");
    i_aq_ = state_variable_indices.at("AEROSOL.A_AQ");
    gas_gas_flat_ = jacobian.VectorIndex(0, i_gas_, i_gas_);
    aq_gas_flat_ = jacobian.VectorIndex(0, i_aq_, i_gas_);
  }

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&,
      DenseMatrixPolicy&) const
  {
  }

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index aq = i_aq_;
    const micm::Real k = rate_constant_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [k](micm::Real& f_gas, micm::Real& f_aq, const micm::Real& gas_val)
              {
                micm::Real rate = k * gas_val;
                f_gas -= rate;
                f_aq += rate;
              },
              forcing_view.GetColumnView(gas),
              forcing_view.GetColumnView(aq),
              state_view.GetConstColumnView(gas));
        },
        forcing,
        state_variables)(forcing, state_variables);
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractJacobianTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& /*state_variables*/,
      SparseMatrixPolicy& jacobian) const
  {
    const micm::Real k = rate_constant_;
    const micm::Index gg = gas_gas_flat_;
    const micm::Index ag = aq_gas_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
          jacobian_view.ForEachBlock(
              [k](micm::Real& j_gg, micm::Real& j_ag)
              {
                j_gg -= -k;
                j_ag -= k;
              },
              jacobian_view.GetBlockView(gg),
              jacobian_view.GetBlockView(ag));
        },
        jacobian)(jacobian);
  }

  // Constraint definition

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    if (!total_mass_.has_value())
    {
      return {};
    }
    return { "AEROSOL.A_AQ" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    if (!total_mass_.has_value())
    {
      return {};
    }
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    if (!total_mass_.has_value())
    {
      return {};
    }
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    return { { i_aq, i_gas }, { i_aq, i_aq } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<class SparseMatrixPolicy>
  void FinalizeConstraintSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    if (!total_mass_.has_value())
    {
      return;
    }
    i_gas_ = state_variable_indices.at("A_GAS");
    i_aq_ = state_variable_indices.at("AEROSOL.A_AQ");
    aq_gas_flat_ = jacobian.VectorIndex(0, i_aq_, i_gas_);
    aq_aq_flat_ = jacobian.VectorIndex(0, i_aq_, i_aq_);
  }

  template<class DenseMatrixPolicy>
  void UpdateConstraintStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&,
      DenseMatrixPolicy&) const
  {
  }

  /// G(y) = [A_GAS] + [A_AQ] - total = 0
  template<class DenseMatrixPolicy>
  void AddConstraintResidual(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index aq = i_aq_;
    const micm::Real total = total_mass_.value();
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [total](micm::Real& f_aq, const micm::Real& gas_val, const micm::Real& aq_val)
              { f_aq = gas_val + aq_val - total; },
              forcing_view.GetColumnView(aq),
              state_view.GetConstColumnView(gas),
              state_view.GetConstColumnView(aq));
        },
        forcing,
        state_variables)(forcing, state_variables);
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractConstraintJacobian(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& /*state_variables*/,
      SparseMatrixPolicy& jacobian) const
  {
    const micm::Index ag = aq_gas_flat_;
    const micm::Index aa = aq_aq_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
          jacobian_view.ForEachBlock(
              [](micm::Real& j_ag, micm::Real& j_aa)
              {
                j_ag -= 1.0;
                j_aa -= 1.0;
              },
              jacobian_view.GetBlockView(ag),
              jacobian_view.GetBlockView(aa));
        },
        jacobian)(jacobian);
  }

 private:
  micm::Real rate_constant_;
  std::optional<micm::Real> total_mass_;
  int i_gas_ = -1;
  int i_aq_ = -1;
  int gas_gas_flat_ = -1;
  int aq_gas_flat_ = -1;
  int aq_aq_flat_ = -1;
};

/// A variant of `StubAerosolWithConstraints` that adds solvent species `S`.
///
/// Process: A_GAS → A_AQ, rate = k * [A_GAS] * [S]
/// Process Jacobian row for A_AQ includes (A_AQ, S), but the constraint only
/// declares (A_AQ, A_GAS) and (A_AQ, A_AQ). Exercises the case where an external
/// model process has a Jacobian element in an algebraic row not declared by any
/// constraint element.
class StubAerosolWithSolvent
{
 public:
  StubAerosolWithSolvent() = delete;

  StubAerosolWithSolvent(micm::Real rate_constant, micm::Real total_mass)
      : rate_constant_(rate_constant),
        total_mass_(total_mass)
  {
  }

  // State definition

  std::tuple<micm::Index, micm::Index> StateSize() const
  {
    return { 1, 0 };
  }
  std::set<std::string> StateVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }
  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  // Process definition

  std::set<std::string> SpeciesUsed() const
  {
    return { "A_GAS", "AEROSOL.A_AQ", "S" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    std::set<std::pair<micm::Index, micm::Index>> elements;
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    auto i_s = state_indices.at("S");
    elements.insert({ i_gas, i_gas });
    elements.insert({ i_gas, i_s });
    elements.insert({ i_aq, i_gas });
    elements.insert({ i_aq, i_s });
    return elements;
  }

  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    i_gas_ = state_variable_indices.at("A_GAS");
    i_aq_ = state_variable_indices.at("AEROSOL.A_AQ");
    i_s_ = state_variable_indices.at("S");
    gas_gas_flat_ = jacobian.VectorIndex(0, i_gas_, i_gas_);
    gas_s_flat_ = jacobian.VectorIndex(0, i_gas_, i_s_);
    aq_gas_flat_ = jacobian.VectorIndex(0, i_aq_, i_gas_);
    aq_s_flat_ = jacobian.VectorIndex(0, i_aq_, i_s_);
  }

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&,
      DenseMatrixPolicy&) const
  {
  }

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index aq = i_aq_;
    const micm::Index s = i_s_;
    const micm::Real k = rate_constant_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [k](
                  micm::Real& f_gas,
                  micm::Real& f_aq,
                  const micm::Real& gas_val,
                  const micm::Real& s_val)
              {
                micm::Real rate = k * gas_val * s_val;
                f_gas -= rate;
                f_aq += rate;
              },
              forcing_view.GetColumnView(gas),
              forcing_view.GetColumnView(aq),
              state_view.GetConstColumnView(gas),
              state_view.GetConstColumnView(s));
        },
        forcing,
        state_variables)(forcing, state_variables);
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractJacobianTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index s = i_s_;
    const micm::Real k = rate_constant_;
    const micm::Index gg = gas_gas_flat_;
    const micm::Index gs = gas_s_flat_;
    const micm::Index ag = aq_gas_flat_;
    const micm::Index as = aq_s_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename SparseMatrixPolicy::ViewType& jacobian_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          jacobian_view.ForEachBlock(
              [k](
                  micm::Real& j_gg,
                  micm::Real& j_gs,
                  micm::Real& j_ag,
                  micm::Real& j_as,
                  const micm::Real& gas_val,
                  const micm::Real& s_val)
              {
                j_gg -= -k * s_val;
                j_ag -= k * s_val;
                j_gs -= -k * gas_val;
                j_as -= k * gas_val;
              },
              jacobian_view.GetBlockView(gg),
              jacobian_view.GetBlockView(gs),
              jacobian_view.GetBlockView(ag),
              jacobian_view.GetBlockView(as),
              state_view.GetConstColumnView(gas),
              state_view.GetConstColumnView(s));
        },
        jacobian,
        state_variables)(jacobian, state_variables);
  }

  // Constraint definition

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }
  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    return { { i_aq, i_gas }, { i_aq, i_aq } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<class SparseMatrixPolicy>
  void FinalizeConstraintSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    i_gas_ = state_variable_indices.at("A_GAS");
    i_aq_ = state_variable_indices.at("AEROSOL.A_AQ");
    aq_gas_flat_ = jacobian.VectorIndex(0, i_aq_, i_gas_);
    aq_aq_flat_ = jacobian.VectorIndex(0, i_aq_, i_aq_);
  }

  template<class DenseMatrixPolicy>
  void UpdateConstraintStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&,
      DenseMatrixPolicy&) const
  {
  }

  template<class DenseMatrixPolicy>
  void AddConstraintResidual(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index aq = i_aq_;
    const micm::Real total = total_mass_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [total](micm::Real& f_aq, const micm::Real& gas_val, const micm::Real& aq_val)
              { f_aq = gas_val + aq_val - total; },
              forcing_view.GetColumnView(aq),
              state_view.GetConstColumnView(gas),
              state_view.GetConstColumnView(aq));
        },
        forcing,
        state_variables)(forcing, state_variables);
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractConstraintJacobian(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& /*state_variables*/,
      SparseMatrixPolicy& jacobian) const
  {
    const micm::Index ag = aq_gas_flat_;
    const micm::Index aa = aq_aq_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
          jacobian_view.ForEachBlock(
              [](micm::Real& j_ag, micm::Real& j_aa)
              {
                j_ag -= 1.0;
                j_aa -= 1.0;
              },
              jacobian_view.GetBlockView(ag),
              jacobian_view.GetBlockView(aa));
        },
        jacobian)(jacobian);
  }

 private:
  micm::Real rate_constant_;
  micm::Real total_mass_;
  int i_gas_ = -1;
  int i_aq_ = -1;
  int i_s_ = -1;
  int gas_gas_flat_ = -1;
  int gas_s_flat_ = -1;
  int aq_gas_flat_ = -1;
  int aq_s_flat_ = -1;
  int aq_aq_flat_ = -1;
};
