// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

/// @file test/integration/stub_aerosol_with_constraints.hpp
/// @brief Stub aerosol model with constraints for integration testing
#pragma once

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <functional>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

/// A stub external model that provides both processes and constraints.
///
/// Simulates gas-to-aerosol partitioning: A_GAS -> A_AQ with rate k,
/// and optionally enforces mass conservation: [A_GAS] + [A_AQ] = total.
///
/// When `total_mass` is provided, the constraint is active (DAE mode).
/// When `total_mass` is std::nullopt, no constraint is active (ODE mode).
class StubAerosolWithConstraints
{
 public:
  StubAerosolWithConstraints() = delete;

  StubAerosolWithConstraints(double rate_constant, std::optional<double> total_mass = std::nullopt)
      : rate_constant_(rate_constant),
        total_mass_(total_mass)
  {
  }

  // ──── State definition methods ────

  std::tuple<std::size_t, std::size_t> StateSize() const
  {
    return { 1, 0 };  // 1 state variable (A_AQ), 0 parameters
  }

  std::set<std::string> StateVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }

  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  // ──── Process definition methods ────

  std::set<std::string> SpeciesUsed() const
  {
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    auto i_gas = state_indices.find("A_GAS");
    auto i_aq = state_indices.find("AEROSOL.A_AQ");
    if (i_gas != state_indices.end() && i_aq != state_indices.end())
    {
      // d(A_GAS)/d(A_GAS) and d(A_AQ)/d(A_GAS)
      elements.insert({ i_gas->second, i_gas->second });
      elements.insert({ i_aq->second, i_gas->second });
    }
    return elements;
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    double k = rate_constant_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        double rate = k * state[i][i_gas];
        forcing[i][i_gas] -= rate;
        forcing[i][i_aq] += rate;
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const SparseMatrixPolicy& jacobian) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    double k = rate_constant_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i_block = 0; i_block < jac.NumberOfBlocks(); ++i_block)
      {
        jac[i_block][i_gas][i_gas] -= (-k);
        jac[i_block][i_aq][i_gas] -= k;
      }
    };
  }

  // ──── Constraint definition methods (optional — satisfies HasConstraints concept) ────

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    if (!total_mass_.has_value())
      return {};
    return { "AEROSOL.A_AQ" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    if (!total_mass_.has_value())
      return {};
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    if (!total_mass_.has_value())
      return {};
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    // Constraint row is i_aq, depends on both A_GAS and A_AQ
    return { { i_aq, i_gas }, { i_aq, i_aq } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  /// Constraint residual: G(y) = [A_GAS] + [A_AQ] - total = 0
  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    double total = total_mass_.value();
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        forcing[i][i_aq] = state[i][i_gas] + state[i][i_aq] - total;
      }
    };
  }

  /// Constraint Jacobian: dG/d[A_GAS] = 1, dG/d[A_AQ] = 1
  /// Subtracted per solver convention: jacobian -= dG/dy
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const SparseMatrixPolicy& jacobian) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i_block = 0; i_block < jac.NumberOfBlocks(); ++i_block)
      {
        jac[i_block][i_aq][i_gas] -= 1.0;
        jac[i_block][i_aq][i_aq] -= 1.0;
      }
    };
  }

 private:
  double rate_constant_;
  std::optional<double> total_mass_;
};

/// A variant of StubAerosolWithConstraints that adds a solvent species `S` dependency.
///
/// Process: A_GAS → A_AQ, rate = k * [A_GAS] * [S]
/// Process Jacobian row for A_AQ includes (A_AQ, S), but the constraint only
/// declares (A_AQ, A_GAS) and (A_AQ, A_AQ).
///
/// This exercises the case where an external model process has a Jacobian element
/// in an algebraic row that is NOT also declared by a constraint element.
class StubAerosolWithSolvent
{
 public:
  StubAerosolWithSolvent() = delete;

  StubAerosolWithSolvent(double rate_constant, double total_mass)
      : rate_constant_(rate_constant),
        total_mass_(total_mass)
  {
  }

  // ──── State definition methods ────

  std::tuple<std::size_t, std::size_t> StateSize() const
  {
    return { 1, 0 };  // 1 state variable (A_AQ), 0 parameters
  }

  std::set<std::string> StateVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }

  std::set<std::string> StateParameterNames() const
  {
    return {};
  }

  // ──── Process definition methods ────

  std::set<std::string> SpeciesUsed() const
  {
    return { "A_GAS", "AEROSOL.A_AQ", "S" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    auto i_s = state_indices.at("S");
    // d(A_GAS)/d(A_GAS), d(A_GAS)/d(S)
    elements.insert({ i_gas, i_gas });
    elements.insert({ i_gas, i_s });
    // d(A_AQ)/d(A_GAS), d(A_AQ)/d(S)  <-- (A_AQ, S) is NOT in constraint elements
    elements.insert({ i_aq, i_gas });
    elements.insert({ i_aq, i_s });
    return elements;
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>& state_parameter_indices) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    auto i_s = state_variable_indices.at("S");
    double k = rate_constant_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        double rate = k * state[i][i_gas] * state[i][i_s];
        forcing[i][i_gas] -= rate;
        forcing[i][i_aq] += rate;
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const SparseMatrixPolicy& jacobian) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    auto i_s = state_variable_indices.at("S");
    double k = rate_constant_;
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jac)
    {
      for (std::size_t i_block = 0; i_block < jac.NumberOfBlocks(); ++i_block)
      {
        double s_val = state[i_block][i_s];
        double gas_val = state[i_block][i_gas];
        // d(rate)/d(A_GAS) = k * [S]
        jac[i_block][i_gas][i_gas] -= (-k * s_val);
        jac[i_block][i_aq][i_gas] -= k * s_val;
        // d(rate)/d(S) = k * [A_GAS]
        jac[i_block][i_gas][i_s] -= (-k * gas_val);
        jac[i_block][i_aq][i_s] -= k * gas_val;   // accesses (A_AQ, S) — triggers bug without fix
      }
    };
  }

  // ──── Constraint definition methods ────

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { "AEROSOL.A_AQ" };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { "A_GAS", "AEROSOL.A_AQ" };
  }

  std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, std::size_t>& state_indices) const
  {
    auto i_gas = state_indices.at("A_GAS");
    auto i_aq = state_indices.at("AEROSOL.A_AQ");
    // Constraint row is i_aq, depends on A_GAS and A_AQ only — NOT S
    return { { i_aq, i_gas }, { i_aq, i_aq } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
      const std::unordered_map<std::string, std::size_t>&) const
  {
    return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    double total = total_mass_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (std::size_t i = 0; i < state.NumRows(); ++i)
      {
        forcing[i][i_aq] = state[i][i_gas] + state[i][i_aq] - total;
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, std::size_t>&,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const SparseMatrixPolicy& jacobian) const
  {
    auto i_gas = state_variable_indices.at("A_GAS");
    auto i_aq = state_variable_indices.at("AEROSOL.A_AQ");
    return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (std::size_t i_block = 0; i_block < jac.NumberOfBlocks(); ++i_block)
      {
        jac[i_block][i_aq][i_gas] -= 1.0;
        jac[i_block][i_aq][i_aq] -= 1.0;
      }
    };
  }

 private:
  double rate_constant_;
  double total_mass_;
};
