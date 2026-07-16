// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Constraint-only external model implementing the Robertson QSSA on species B:
//   G(A,B,C) = k1*A - k2*B^2 - k3*B*C = 0   (replaces dB/dt = 0)
//
// B, A, C are pre-existing species, so this model declares no new state
// variables and is wired via SolverBuilder::AddExternalModel only (no entry in
// System's external-model collection).
#pragma once

#include <micm/system/conditions.hpp>

#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace robertson
{
  class QssaConstraint
  {
   public:
    QssaConstraint(double k1, double k2, double k3)
        : k1_(k1),
          k2_(k2),
          k3_(k3)
    {
    }

    // B's ODE row becomes algebraic.
    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "B" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "A", "B", "C" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& state_indices) const
    {
      auto i_a = state_indices.at("A");
      auto i_b = state_indices.at("B");
      auto i_c = state_indices.at("C");
      // Constraint row is B; depends on A, B, C.
      return { { i_b, i_a }, { i_b, i_b }, { i_b, i_c } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>
    ConstraintUpdateStateParametersFunction(const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
    }

    // Residual: G = k1*A - k2*B^2 - k3*B*C  (written into B's forcing row).
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& var) const
    {
      auto i_a = var.at("A");
      auto i_b = var.at("B");
      auto i_c = var.at("C");
      double k1 = k1_, k2 = k2_, k3 = k3_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t i = 0; i < state.NumRows(); ++i)
        {
          double a = state[i][i_a], b = state[i][i_b], c = state[i][i_c];
          forcing[i][i_b] = k1 * a - k2 * b * b - k3 * b * c;
        }
      };
    }

    // Jacobian (solver subtracts dG/dy):
    //   dG/dA = k1 ; dG/dB = -2*k2*B - k3*C ; dG/dC = -k3*B
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& var,
        const SparseMatrixPolicy&) const
    {
      auto i_a = var.at("A");
      auto i_b = var.at("B");
      auto i_c = var.at("C");
      double k1 = k1_, k2 = k2_, k3 = k3_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
      {
        for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
        {
          double b = state[i][i_b], c = state[i][i_c];
          jac[i][i_b][i_a] -= k1;
          jac[i][i_b][i_b] -= (-2.0 * k2 * b - k3 * c);
          jac[i][i_b][i_c] -= (-k3 * b);
        }
      };
    }

   private:
    double k1_;
    double k2_;
    double k3_;
  };
}  // namespace robertson
