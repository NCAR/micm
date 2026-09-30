// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_cloud_chemistry.hpp
/// @brief Stub cloud chemistry model for integration testing
#pragma once

#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <cmath>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

/// A stub external model with gas-aqueous partitioning and aqueous chemistry, like a cloud model.
///
/// All concentrations are in mol m-3 of air.
///
/// Constraint (A_AQ row is algebraic): Henry's law equilibrium
///   G = K_p(T) * [A_G] - [A_AQ] = 0
///   K_p(T) = H(T) * R * T * L    (dimensionless partition coefficient)
///   H(T) = H_ref * exp(dH_R * (1/T - 1/T_ref))    (van 't Hoff, mol L-1 atm-1)
///   L = LWC / rho_water    (liquid water volume fraction, m3 water per m3 air)
///
/// Process: aqueous reaction A_AQ -> P_AQ at rate k * [A_AQ]
///
/// The process writes into the A_AQ row, which the constraint owns. This is the structure
/// that https://github.com/NCAR/micm/issues/1094 addresses.
class StubCloudChemistry
{
 public:
  /// Gas constant in L atm mol-1 K-1, to match the units of the Henry's law constant
  static constexpr micm::Real GAS_CONSTANT = 0.082057;
  /// Density of liquid water in kg m-3
  static constexpr micm::Real WATER_DENSITY = 1000.0;

  StubCloudChemistry() = delete;

  /// @param rate_constant Aqueous rate constant k (s-1)
  /// @param henry_constant_ref Henry's law constant at T_ref (mol L-1 atm-1)
  /// @param delta_h_over_r Temperature dependence of the Henry's law constant (K)
  /// @param liquid_water_content Liquid water content (kg m-3)
  /// @param temperature_ref Reference temperature (K)
  StubCloudChemistry(
      micm::Real rate_constant,
      micm::Real henry_constant_ref,
      micm::Real delta_h_over_r,
      micm::Real liquid_water_content,
      micm::Real temperature_ref = 298.15)
      : rate_constant_(rate_constant),
        henry_constant_ref_(henry_constant_ref),
        delta_h_over_r_(delta_h_over_r),
        liquid_water_content_(liquid_water_content),
        temperature_ref_(temperature_ref)
  {
  }

  /// @brief Returns the partition coefficient K_p at a temperature. The tests use it for the exact solution.
  micm::Real PartitionCoefficient(micm::Real temperature) const
  {
    return PartitionCoefficient(
        temperature, henry_constant_ref_, delta_h_over_r_, temperature_ref_, liquid_water_content_ / WATER_DENSITY);
  }

  // State definition

  std::tuple<micm::Index, micm::Index> StateSize() const
  {
    return { 2, 0 };
  }
  std::set<std::string> StateVariableNames() const
  {
    return { "CLOUD.A_AQ", "CLOUD.P_AQ" };
  }
  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  // Process definition

  std::set<std::string> SpeciesUsed() const
  {
    return { "CLOUD.A_AQ", "CLOUD.P_AQ" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    auto i_aq = state_indices.at("CLOUD.A_AQ");
    auto i_p = state_indices.at("CLOUD.P_AQ");
    return { { i_aq, i_aq }, { i_p, i_aq } };
  }

  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& /*state_parameter_indices*/,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    i_gas_ = state_variable_indices.at("A_G");
    i_aq_ = state_variable_indices.at("CLOUD.A_AQ");
    i_p_ = state_variable_indices.at("CLOUD.P_AQ");
    aq_aq_flat_ = jacobian.VectorIndex(0, i_aq_, i_aq_);
    p_aq_flat_ = jacobian.VectorIndex(0, i_p_, i_aq_);
  }

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)
      const
  {
  }

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& /*state_parameters*/,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index aq = i_aq_;
    const micm::Index p = i_p_;
    const micm::Real k = rate_constant_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [k](micm::Real& f_aq, micm::Real& f_p, const micm::Real& aq_val)
              {
                const micm::Real rate = k * aq_val;
                f_aq -= rate;
                f_p += rate;
              },
              forcing_view.GetColumnView(aq),
              forcing_view.GetColumnView(p),
              state_view.GetConstColumnView(aq));
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
    const micm::Index aa = aq_aq_flat_;
    const micm::Index pa = p_aq_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(const typename SparseMatrixPolicy::ViewType& jacobian_view) {
          jacobian_view.ForEachBlock(
              [k](micm::Real& j_aa, micm::Real& j_pa)
              {
                j_aa -= -k;
                j_pa -= k;
              },
              jacobian_view.GetBlockView(aa),
              jacobian_view.GetBlockView(pa));
        },
        jacobian)(jacobian);
  }

  // Constraint definition

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { "CLOUD.A_AQ" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { "A_G", "CLOUD.A_AQ" };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    auto i_gas = state_indices.at("A_G");
    auto i_aq = state_indices.at("CLOUD.A_AQ");
    return { { i_aq, i_gas }, { i_aq, i_aq } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return { "CLOUD.K_P" };
  }

  template<class SparseMatrixPolicy>
  void FinalizeConstraintSetup(
      const std::unordered_map<std::string, micm::Index>& state_parameter_indices,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian)
  {
    i_gas_ = state_variable_indices.at("A_G");
    i_aq_ = state_variable_indices.at("CLOUD.A_AQ");
    i_k_p_ = state_parameter_indices.at("CLOUD.K_P");
    aq_gas_flat_ = jacobian.VectorIndex(0, i_aq_, i_gas_);
    aq_aq_flat_ = jacobian.VectorIndex(0, i_aq_, i_aq_);
  }

  /// Computes K_p(T) for each grid cell
  template<class DenseMatrixPolicy>
  void UpdateConstraintStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
      DenseMatrixPolicy& params) const
  {
    const micm::Index i_k_p = i_k_p_;
    const micm::Real h_ref = henry_constant_ref_;
    const micm::Real dh_r = delta_h_over_r_;
    const micm::Real t_ref = temperature_ref_;
    const micm::Real water_fraction = liquid_water_content_ / WATER_DENSITY;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
            const typename DenseMatrixPolicy::ViewType& params_view) {
          params_view.ForEachRow(
              [h_ref, dh_r, t_ref, water_fraction](const micm::Conditions& cond, micm::Real& k_p)
              { k_p = PartitionCoefficient(cond.temperature_, h_ref, dh_r, t_ref, water_fraction); },
              conditions_view,
              params_view.GetColumnView(i_k_p));
        },
        conditions,
        params)(conditions, params);
  }

  /// G = K_p * [A_G] - [A_AQ]
  template<class DenseMatrixPolicy>
  void AddConstraintResidual(
      const DenseMatrixPolicy& params,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const
  {
    const micm::Index gas = i_gas_;
    const micm::Index aq = i_aq_;
    const micm::Index i_k_p = i_k_p_;
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& params_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view) {
          forcing_view.ForEachRow(
              [](micm::Real& f_aq, const micm::Real& k_p, const micm::Real& gas_val, const micm::Real& aq_val)
              { f_aq = k_p * gas_val - aq_val; },
              forcing_view.GetColumnView(aq),
              params_view.GetConstColumnView(i_k_p),
              state_view.GetConstColumnView(gas),
              state_view.GetConstColumnView(aq));
        },
        forcing,
        params,
        state_variables)(forcing, params, state_variables);
  }

  /// Subtracts dG/dy from the Jacobian (solver convention)
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractConstraintJacobian(
      const DenseMatrixPolicy& params,
      const DenseMatrixPolicy& /*state_variables*/,
      SparseMatrixPolicy& jacobian) const
  {
    const micm::Index i_k_p = i_k_p_;
    const micm::Index ag = aq_gas_flat_;
    const micm::Index aa = aq_aq_flat_;
    SparseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename SparseMatrixPolicy::ViewType& jacobian_view,
            const typename DenseMatrixPolicy::ConstViewType& params_view) {
          jacobian_view.ForEachBlock(
              [](micm::Real& j_ag, micm::Real& j_aa, const micm::Real& k_p)
              {
                j_ag -= k_p;
                j_aa -= -1.0;
              },
              jacobian_view.GetBlockView(ag),
              jacobian_view.GetBlockView(aa),
              params_view.GetConstColumnView(i_k_p));
        },
        jacobian,
        params)(jacobian, params);
  }

 private:
  static MICM_INLINE_DEVICE_FUNCTION micm::Real PartitionCoefficient(
      micm::Real temperature,
      micm::Real henry_constant_ref,
      micm::Real delta_h_over_r,
      micm::Real temperature_ref,
      micm::Real water_fraction)
  {
    const micm::Real henry_constant =
        henry_constant_ref * std::exp(delta_h_over_r * (1.0 / temperature - 1.0 / temperature_ref));
    return henry_constant * GAS_CONSTANT * temperature * water_fraction;
  }

  micm::Real rate_constant_;
  micm::Real henry_constant_ref_;
  micm::Real delta_h_over_r_;
  micm::Real liquid_water_content_;
  micm::Real temperature_ref_;
  int i_gas_ = -1;
  int i_aq_ = -1;
  int i_p_ = -1;
  int i_k_p_ = -1;
  int aq_aq_flat_ = -1;
  int p_aq_flat_ = -1;
  int aq_gas_flat_ = -1;
};
