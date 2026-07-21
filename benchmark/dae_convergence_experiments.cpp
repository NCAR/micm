// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Numerical experiments for DAE constraint initialization, convergence order,
// tolerance coupling, nonlinear robustness, and adaptive-controller behavior.

#include <micm/CPU.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  constexpr double TEMPERATURE = 298.0;
  constexpr double PRESSURE = 101325.0;

  void SetConditions(auto& state)
  {
    for (auto& conditions : state.conditions_)
    {
      conditions.temperature_ = TEMPERATURE;
      conditions.pressure_ = PRESSURE;
    }
  }

  std::string StateName(micm::SolverState state)
  {
    return micm::SolverStateToString(state);
  }

  double RelativeError(double value, double reference)
  {
    return std::abs(value - reference) / std::max(std::abs(reference), 1.0e-300);
  }

  std::vector<double> PowersOfTen(int first, int last, int stride)
  {
    std::vector<double> values;
    for (int exponent = first; exponent <= last; exponent += stride)
    {
      values.push_back(std::pow(10.0, exponent));
    }
    return values;
  }

  enum class ManufacturedManifold
  {
    Quadratic,
    Sine
  };

  /// Manufactured index-1 DAE X'=-Z with either a quadratic or sine constraint manifold.
  class ManufacturedDaeModel
  {
   public:
    ManufacturedDaeModel(
        double row_scale,
        double state_scale,
        ManufacturedManifold manifold = ManufacturedManifold::Quadratic)
        : row_scale_(row_scale),
          state_scale_(state_scale),
          manifold_(manifold)
    {
    }

    std::set<std::string> SpeciesUsed() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      return { { indices.at("X"), indices.at("Z") } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][x] -= state[cell][z];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][x][z] -= -1.0;
        }
      };
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "Z" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto z = indices.at("Z");
      return { { z, indices.at("X") }, { z, z } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double row_scale = row_scale_;
      const double state_scale = state_scale_;
      const auto manifold = manifold_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double constrained_value = manifold == ManufacturedManifold::Quadratic
                                               ? state[cell][x] * state[cell][x] / state_scale
                                               : state_scale * std::sin(state[cell][x] / state_scale);
          forcing[cell][z] = row_scale * (state[cell][z] - constrained_value);
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double row_scale = row_scale_;
      const double state_scale = state_scale_;
      const auto manifold = manifold_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double derivative = manifold == ManufacturedManifold::Quadratic ? 2.0 * state[block][x] / state_scale
                                                                                : std::cos(state[block][x] / state_scale);
          jacobian[block][z][x] -= -row_scale * derivative;
          jacobian[block][z][z] -= row_scale;
        }
      };
    }

   private:
    double row_scale_;
    double state_scale_;
    ManufacturedManifold manifold_;
  };

  /// Constraint-only model G = row_scale * (Z^2 - X).
  class SquareRootConstraintModel
  {
   public:
    explicit SquareRootConstraintModel(double row_scale)
        : row_scale_(row_scale)
    {
    }

    std::set<std::string> SpeciesUsed() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      return { { indices.at("X"), indices.at("Z") } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][x] -= state[cell][z];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][x][z] -= -1.0;
        }
      };
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "Z" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto z = indices.at("Z");
      return { { z, indices.at("X") }, { z, z } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double row_scale = row_scale_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][z] = row_scale * (state[cell][z] * state[cell][z] - state[cell][x]);
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double row_scale = row_scale_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][z][x] -= -row_scale;
          jacobian[block][z][z] -= 2.0 * row_scale * state[block][z];
        }
      };
    }

   private:
    double row_scale_;
  };

  /// Two linear constraints with a tunable nearly singular algebraic block.
  class IllConditionedConstraintModel
  {
   public:
    explicit IllConditionedConstraintModel(double epsilon)
        : epsilon_(epsilon)
    {
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "Z1", "Z2" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "Z1", "Z2" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto z1 = indices.at("Z1");
      const auto z2 = indices.at("Z2");
      return { { z1, z1 }, { z1, z2 }, { z2, z1 }, { z2, z2 } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto z1 = indices.at("Z1");
      const auto z2 = indices.at("Z2");
      const double epsilon = epsilon_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][z1] = state[cell][z1] + state[cell][z2] - 2.0;
          forcing[cell][z2] = state[cell][z1] + (1.0 + epsilon) * state[cell][z2] - (2.0 + epsilon);
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto z1 = indices.at("Z1");
      const auto z2 = indices.at("Z2");
      const double epsilon = epsilon_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][z1][z1] -= 1.0;
          jacobian[block][z1][z2] -= 1.0;
          jacobian[block][z2][z1] -= 1.0;
          jacobian[block][z2][z2] -= 1.0 + epsilon;
        }
      };
    }

   private:
    double epsilon_;
  };

  /// X'=-X with any number of algebraic copies Z_j=X.
  class SlavedVariablesModel
  {
   public:
    explicit SlavedVariablesModel(std::size_t number_of_algebraic_variables, bool constant_forcing = false)
        : constant_forcing_(constant_forcing)
    {
      for (std::size_t i = 0; i < number_of_algebraic_variables; ++i)
      {
        algebraic_names_.push_back("Z" + std::to_string(i));
      }
    }

    std::set<std::string> SpeciesUsed() const
    {
      std::set<std::string> names{ "X" };
      names.insert(algebraic_names_.begin(), algebraic_names_.end());
      return names;
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      if (constant_forcing_)
      {
        return {};
      }
      const auto x = indices.at("X");
      return { { x, x } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const bool constant_forcing = constant_forcing_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][x] -= constant_forcing ? 1.0 : state[cell][x];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const bool constant_forcing = constant_forcing_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        if (constant_forcing)
        {
          return;
        }
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][x][x] -= -1.0;
        }
      };
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { algebraic_names_.begin(), algebraic_names_.end() };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      std::set<std::string> names{ "X" };
      names.insert(algebraic_names_.begin(), algebraic_names_.end());
      return names;
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      std::set<std::pair<std::size_t, std::size_t>> elements;
      for (const auto& name : algebraic_names_)
      {
        const auto z = indices.at(name);
        elements.insert({ z, x });
        elements.insert({ z, z });
      }
      return elements;
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      std::vector<std::size_t> algebraic_indices;
      for (const auto& name : algebraic_names_)
      {
        algebraic_indices.push_back(indices.at(name));
      }
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          for (const auto z : algebraic_indices)
          {
            forcing[cell][z] = state[cell][z] - state[cell][x];
          }
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      std::vector<std::size_t> algebraic_indices;
      for (const auto& name : algebraic_names_)
      {
        algebraic_indices.push_back(indices.at(name));
      }
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          for (const auto z : algebraic_indices)
          {
            jacobian[block][z][x] -= -1.0;
            jacobian[block][z][z] -= 1.0;
          }
        }
      };
    }

   private:
    std::vector<std::string> algebraic_names_;
    bool constant_forcing_;
  };

  /// Van der Pol oscillator in singular-perturbation (Lienard) form, with the
  /// fast variable negated (Z := -dX/dt) so the whole trajectory stays in the
  /// positive orthant — MICM clamps species nonnegative on exit, which would
  /// otherwise zero the reported fast variable:
  ///   X' = -Z,   eps * Z' = (1 - X^2) Z + X.
  /// The eps = 0 limit is the index-1 DAE with G = (1 - X^2) Z + X = 0 on the
  /// Z row (Z = X / (X^2 - 1), positive for X > 1) — the same structural limit
  /// QSSA takes for a fast species. The epsilon experiment measures uniform
  /// convergence of the stiff-ODE terminal state to the DAE terminal state as
  /// eps -> 0.
  ///
  /// Two model classes so the eps > 0 case is a genuinely constraint-free
  /// external model (the HasConstraints concept is not satisfied).
  class VanDerPolOdeModel
  {
   public:
    explicit VanDerPolOdeModel(double epsilon)
        : epsilon_(epsilon)
    {
    }

    std::set<std::string> SpeciesUsed() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return { { x, z }, { z, x }, { z, z } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double epsilon = epsilon_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double xv = state[cell][x];
          const double zv = state[cell][z];
          forcing[cell][x] -= zv;
          forcing[cell][z] += ((1.0 - xv * xv) * zv + xv) / epsilon;
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      const double epsilon = epsilon_;
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double xv = state[block][x];
          const double zv = state[block][z];
          jacobian[block][x][z] -= -1.0;
          jacobian[block][z][x] -= (-2.0 * xv * zv + 1.0) / epsilon;
          jacobian[block][z][z] -= (1.0 - xv * xv) / epsilon;
        }
      };
    }

   private:
    double epsilon_;
  };

  /// The eps = 0 limit: X row differential (X' = -Z), Z row algebraic
  /// (G = (1 - X^2) Z + X = 0).
  class VanDerPolDaeModel
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      return { { indices.at("X"), indices.at("Z") } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][x] -= state[cell][z];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][x][z] -= -1.0;
        }
      };
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "Z" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "X", "Z" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto z = indices.at("Z");
      return { { z, indices.at("X") }, { z, z } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double xv = state[cell][x];
          forcing[cell][z] = (1.0 - xv * xv) * state[cell][z] + xv;
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& indices,
        const SparseMatrixPolicy&) const
    {
      const auto x = indices.at("X");
      const auto z = indices.at("Z");
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double xv = state[block][x];
          const double zv = state[block][z];
          jacobian[block][z][x] -= (-2.0 * xv * zv + 1.0);
          jacobian[block][z][z] -= (1.0 - xv * xv);
        }
      };
    }
  };

  auto BuildManufacturedSolver(
      const micm::RosenbrockSolverParameters& parameters,
      double row_scale,
      double state_scale,
      ManufacturedManifold manifold = ManufacturedManifold::Quadratic)
  {
    const auto x = micm::Species("X");
    const auto z = micm::Species("Z");
    const micm::Phase gas{ "gas", std::vector<micm::PhaseSpecies>{ x, z } };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(ManufacturedDaeModel(row_scale, state_scale, manifold))
        .SetReorderState(false)
        .Build();
  }

  auto BuildSquareRootSolver(const micm::RosenbrockSolverParameters& parameters, double row_scale)
  {
    const auto x = micm::Species("X");
    const auto z = micm::Species("Z");
    const micm::Phase gas{ "gas", std::vector<micm::PhaseSpecies>{ x, z } };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(SquareRootConstraintModel(row_scale))
        .SetReorderState(false)
        .Build();
  }

  auto BuildIllConditionedSolver(const micm::RosenbrockSolverParameters& parameters, double epsilon)
  {
    const auto z1 = micm::Species("Z1");
    const auto z2 = micm::Species("Z2");
    const micm::Phase gas{ "gas", std::vector<micm::PhaseSpecies>{ z1, z2 } };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(IllConditionedConstraintModel(epsilon))
        .SetReorderState(false)
        .Build();
  }

  auto BuildSlavedSolver(
      const micm::RosenbrockSolverParameters& parameters,
      std::size_t algebraic_variables,
      bool constant_forcing = false)
  {
    std::vector<micm::PhaseSpecies> species{ micm::Species("X") };
    for (std::size_t i = 0; i < algebraic_variables; ++i)
    {
      species.push_back(micm::Species("Z" + std::to_string(i)));
    }
    const micm::Phase gas{ "gas", species };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(SlavedVariablesModel(algebraic_variables, constant_forcing))
        .SetReorderState(false)
        .Build();
  }

  auto BuildVanDerPolOdeSolver(const micm::RosenbrockSolverParameters& parameters, double epsilon)
  {
    const auto x = micm::Species("X");
    const auto z = micm::Species("Z");
    const micm::Phase gas{ "gas", std::vector<micm::PhaseSpecies>{ x, z } };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(VanDerPolOdeModel(epsilon))
        .SetReorderState(false)
        .Build();
  }

  auto BuildVanDerPolDaeSolver(const micm::RosenbrockSolverParameters& parameters)
  {
    const auto x = micm::Species("X");
    const auto z = micm::Species("Z");
    const micm::Phase gas{ "gas", std::vector<micm::PhaseSpecies>{ x, z } };
    return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(parameters)
        .SetSystem(micm::System(gas))
        .AddExternalModel(VanDerPolDaeModel())
        .SetReorderState(false)
        .Build();
  }

  void RunProjectionScaling(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "projection_scaling.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "x,initial_relative_error,row_scale,constraint_tolerance,max_iterations,state,iterations,"
           "z_initial,z_final,z_exact,relative_state_error,raw_residual,unscaled_residual,weighted_initial_correction\n";

    const std::array<double, 3> x_values{ 1.0e-6, 1.0, 1.0e6 };
    const std::array<double, 5> initial_errors{ 1.0e-12, 1.0e-6, 1.0e-2, 1.0, 1.0e6 };
    const auto row_scales = PowersOfTen(-14, 14, 2);
    const std::array<double, 6> tolerances{ 1.0, 1.0e-1, 1.0e-2, 1.0e-4, 1.0e-8, 1.0e-12 };
    const std::array<std::size_t, 4> max_iterations{ 0, 1, 2, 10 };

    std::size_t false_convergence = 0;
    std::size_t solved_but_failed = 0;
    std::size_t cases = 0;

    for (const double row_scale : row_scales)
    {
      auto base_parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      auto solver = BuildManufacturedSolver(base_parameters, row_scale, 1.0);
      for (const double x_value : x_values)
      {
        const double z_exact = x_value * x_value;
        for (const double initial_error : initial_errors)
        {
          const double z_initial = z_exact * (1.0 + initial_error);
          for (const double tolerance : tolerances)
          {
            for (const std::size_t maximum : max_iterations)
            {
              auto parameters = base_parameters;
              parameters.constraint_init_tolerance_ = tolerance;
              parameters.constraint_init_max_iterations_ = maximum;

              auto state = solver.GetState(1);
              SetConditions(state);
              state.SetRelativeTolerance(1.0e-6);
              state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-12));
              const auto x = state.variable_map_.at("X");
              const auto z = state.variable_map_.at("Z");
              state.variables_[0][x] = x_value;
              state.variables_[0][z] = z_initial;

              micm::SolverStats stats;
              const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
              const double z_final = state.variables_[0][z];
              const double relative_state_error = RelativeError(z_final, z_exact);
              const double unscaled_residual = z_final - z_exact;
              const double raw_residual = row_scale * unscaled_residual;
              const double initial_correction = std::abs(z_initial - z_exact);
              const double correction_scale = 1.0e-12 + 1.0e-6 * std::max(std::abs(z_initial), std::abs(z_exact));
              const double weighted_initial_correction = initial_correction / correction_scale;

              if (status == micm::SolverState::Converged && relative_state_error > 1.0e-6)
              {
                ++false_convergence;
              }
              if (status != micm::SolverState::Converged && relative_state_error <= 1.0e-12)
              {
                ++solved_but_failed;
              }
              ++cases;

              csv << x_value << ',' << initial_error << ',' << row_scale << ',' << tolerance << ',' << maximum << ','
                  << StateName(status) << ',' << stats.constraint_init_iterations_ << ',' << z_initial << ',' << z_final
                  << ',' << z_exact << ',' << relative_state_error << ',' << raw_residual << ',' << unscaled_residual << ','
                  << weighted_initial_correction << '\n';
            }
          }
        }
      }
    }

    summary << "- Projection scaling: " << cases << " cases; false convergence=" << false_convergence
            << ", solved-but-reported-failed=" << solved_but_failed << ".\n";
  }

  void RunFixedStepOrder(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "fixed_step_order.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "method,expected_order,k,H,row_scale,state_scale,state,steps,rejected,x_error,z_error,raw_residual,"
           "observed_x_order,observed_z_order\n";

    struct Method
    {
      const char* name;
      double expected_order;
      std::function<micm::RosenbrockSolverParameters()> parameters;
    };
    const std::array<Method, 2> methods{
      Method{ "four_stage_dae",
              3.0,
              []() { return micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters(); } },
      Method{ "six_stage_dae",
              4.0,
              []() { return micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters(); } },
    };
    const std::array<double, 3> row_scales{ 1.0e-12, 1.0, 1.0e12 };
    const std::array<double, 3> state_scales{ 1.0e-12, 1.0, 1.0e12 };

    for (const auto& method : methods)
    {
      for (const double row_scale : row_scales)
      {
        for (const double state_scale : state_scales)
        {
          double previous_x_error = std::numeric_limits<double>::quiet_NaN();
          double previous_z_error = std::numeric_limits<double>::quiet_NaN();
          for (int k = 3; k <= 11; ++k)
          {
            const double H = std::ldexp(1.0, -k);
            auto parameters = method.parameters();
            parameters.h_start_ = H;
            parameters.h_min_ = H;
            parameters.h_max_ = H;
            parameters.max_number_of_steps_ = 10000000;
            auto solver = BuildManufacturedSolver(parameters, row_scale, state_scale, ManufacturedManifold::Sine);
            auto state = solver.GetState(1);
            SetConditions(state);
            state.SetRelativeTolerance(1.0e6);
            state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e6 * state_scale));
            const auto x = state.variable_map_.at("X");
            const auto z = state.variable_map_.at("Z");
            state.variables_[0][x] = state_scale;
            state.variables_[0][z] = state_scale * std::sin(1.0);

            const auto result = solver.Solve(1.0, state, parameters);
            const double x_exact = 2.0 * state_scale * std::atan(std::tan(0.5) * std::exp(-1.0));
            const double z_exact = state_scale * std::sin(x_exact / state_scale);
            const double x_error = RelativeError(state.variables_[0][x], x_exact);
            const double z_error = RelativeError(state.variables_[0][z], z_exact);
            const double raw_residual =
                row_scale * (state.variables_[0][z] - state_scale * std::sin(state.variables_[0][x] / state_scale));
            const double x_order = std::isfinite(previous_x_error) && x_error > 0.0
                                       ? std::log2(previous_x_error / x_error)
                                       : std::numeric_limits<double>::quiet_NaN();
            const double z_order = std::isfinite(previous_z_error) && z_error > 0.0
                                       ? std::log2(previous_z_error / z_error)
                                       : std::numeric_limits<double>::quiet_NaN();

            csv << method.name << ',' << method.expected_order << ',' << k << ',' << H << ',' << row_scale << ','
                << state_scale << ',' << StateName(result.state_) << ',' << result.stats_.accepted_ << ','
                << result.stats_.rejected_ << ',' << x_error << ',' << z_error << ',' << raw_residual << ',' << x_order
                << ',' << z_order << '\n';

            previous_x_error = x_error;
            previous_z_error = z_error;
          }
        }
      }
    }

    summary << "- Fixed-step order: completed both DAE methods over H=2^-3...2^-11 and 3x3 row/state scales.\n";
  }

  void RunToleranceCoupling(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "tolerance_coupling.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "row_scale,rtol,atol,constraint_tolerance,initial_relative_error,state,init_iterations,steps,rejected,"
           "x_relative_error,z_relative_error,raw_residual,unscaled_residual\n";

    const std::array<double, 3> row_scales{ 1.0e-8, 1.0, 1.0e8 };
    const std::array<double, 4> rtols{ 1.0e-3, 1.0e-5, 1.0e-7, 1.0e-9 };
    const std::array<double, 6> constraint_tolerances{ 1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-6 };
    const std::array<double, 4> initial_errors{ 1.0e-6, 1.0e-2, 1.0, 1.0e2 };
    std::size_t failures = 0;

    for (const double row_scale : row_scales)
    {
      auto base_parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      base_parameters.max_number_of_steps_ = 100000;
      base_parameters.constraint_init_max_iterations_ = 20;
      auto solver = BuildSquareRootSolver(base_parameters, row_scale);
      for (const double rtol : rtols)
      {
        const double atol = rtol * 1.0e-3;
        for (const double constraint_tolerance : constraint_tolerances)
        {
          for (const double initial_error : initial_errors)
          {
            auto parameters = base_parameters;
            parameters.constraint_init_tolerance_ = constraint_tolerance;
            auto state = solver.GetState(1);
            SetConditions(state);
            state.SetRelativeTolerance(rtol);
            state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, atol));
            const auto x = state.variable_map_.at("X");
            const auto z = state.variable_map_.at("Z");
            state.variables_[0][x] = 1.0;
            state.variables_[0][z] = 1.0 + initial_error;

            const auto result = solver.Solve(1.0, state, parameters);
            if (result.state_ != micm::SolverState::Converged)
            {
              ++failures;
            }
            const double x_error = RelativeError(state.variables_[0][x], 0.25);
            const double z_error = RelativeError(state.variables_[0][z], 0.5);
            const double unscaled_residual = state.variables_[0][z] * state.variables_[0][z] - state.variables_[0][x];

            csv << row_scale << ',' << rtol << ',' << atol << ',' << constraint_tolerance << ',' << initial_error << ','
                << StateName(result.state_) << ',' << result.stats_.constraint_init_iterations_ << ','
                << result.stats_.accepted_ << ',' << result.stats_.rejected_ << ',' << x_error << ',' << z_error << ','
                << row_scale * unscaled_residual << ',' << unscaled_residual << '\n';
          }
        }
      }
    }

    summary << "- Tolerance coupling: 288 cases; non-converged=" << failures << ".\n";
  }

  void RunNonlinearInitialization(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "nonlinear_initialization.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "problem,row_scale_or_epsilon,initial_z1,initial_z2,state,iterations,final_z1,final_z2,"
           "state_error,residual_1,residual_2\n";

    const std::array<double, 3> row_scales{ 1.0e-8, 1.0, 1.0e8 };
    const std::array<double, 9> initial_z{ -10.0, -1.0, -1.0e-6, 0.0, 1.0e-6, 0.1, 1.0, 2.0, 10.0 };
    std::size_t nonlinear_failures = 0;
    for (const double row_scale : row_scales)
    {
      auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      parameters.constraint_init_tolerance_ = 0.1;
      auto solver = BuildSquareRootSolver(parameters, row_scale);
      for (const double initial : initial_z)
      {
        auto state = solver.GetState(1);
        SetConditions(state);
        const auto x = state.variable_map_.at("X");
        const auto z = state.variable_map_.at("Z");
        state.variables_[0][x] = 1.0;
        state.variables_[0][z] = initial;
        micm::SolverStats stats;
        const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
        const double final_z = state.variables_[0][z];
        const double residual = row_scale * (final_z * final_z - 1.0);
        if (status != micm::SolverState::Converged)
        {
          ++nonlinear_failures;
        }
        csv << "sqrt," << row_scale << ',' << initial << ",nan," << StateName(status) << ','
            << stats.constraint_init_iterations_ << ',' << final_z << ",nan," << std::abs(final_z - 1.0) << ',' << residual
            << ",nan\n";
      }
    }

    std::size_t conditioning_failures = 0;
    for (int exponent = -2; exponent >= -14; exponent -= 2)
    {
      const double epsilon = std::pow(10.0, exponent);
      auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      parameters.constraint_init_tolerance_ = 0.1;
      auto solver = BuildIllConditionedSolver(parameters, epsilon);
      auto state = solver.GetState(1);
      SetConditions(state);
      const auto z1 = state.variable_map_.at("Z1");
      const auto z2 = state.variable_map_.at("Z2");
      state.variables_[0][z1] = 0.0;
      state.variables_[0][z2] = 0.0;
      micm::SolverStats stats;
      const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
      const double final_z1 = state.variables_[0][z1];
      const double final_z2 = state.variables_[0][z2];
      const double residual_1 = final_z1 + final_z2 - 2.0;
      const double residual_2 = final_z1 + (1.0 + epsilon) * final_z2 - (2.0 + epsilon);
      const double error = std::max(std::abs(final_z1 - 1.0), std::abs(final_z2 - 1.0));
      if (status != micm::SolverState::Converged)
      {
        ++conditioning_failures;
      }
      csv << "ill_conditioned," << epsilon << ",0,0," << StateName(status) << ',' << stats.constraint_init_iterations_ << ','
          << final_z1 << ',' << final_z2 << ',' << error << ',' << residual_1 << ',' << residual_2 << '\n';
    }

    summary << "- Nonlinear initialization: sqrt failures=" << nonlinear_failures
            << "; ill-conditioned failures=" << conditioning_failures << ".\n";
  }

  void RunControllerAndCadence(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "controller_cadence.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "experiment,algebraic_variables,cells,solve_calls,rtol,state,steps,rejected,init_iterations,"
           "x_relative_error,max_algebraic_relative_error,max_raw_residual\n";

    const std::array<std::size_t, 4> algebraic_counts{ 0, 1, 9, 99 };
    const std::array<double, 3> rtols{ 1.0e-4, 1.0e-6, 1.0e-8 };
    for (const std::size_t algebraic_count : algebraic_counts)
    {
      auto base_parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      base_parameters.max_number_of_steps_ = 100000;
      auto solver = BuildSlavedSolver(base_parameters, algebraic_count);
      for (const double rtol : rtols)
      {
        auto state = solver.GetState(1);
        SetConditions(state);
        state.SetRelativeTolerance(rtol);
        state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, rtol * 1.0e-3));
        const auto x = state.variable_map_.at("X");
        state.variables_[0][x] = 1.0;
        for (std::size_t i = 0; i < algebraic_count; ++i)
        {
          state.variables_[0][state.variable_map_.at("Z" + std::to_string(i))] = 1.0;
        }
        const auto result = solver.Solve(1.0, state, base_parameters);
        const double exact = std::exp(-1.0);
        double maximum_algebraic_error = 0.0;
        double maximum_residual = 0.0;
        for (std::size_t i = 0; i < algebraic_count; ++i)
        {
          const double value = state.variables_[0][state.variable_map_.at("Z" + std::to_string(i))];
          maximum_algebraic_error = std::max(maximum_algebraic_error, RelativeError(value, exact));
          maximum_residual = std::max(maximum_residual, std::abs(value - state.variables_[0][x]));
        }
        csv << "algebraic_fraction," << algebraic_count << ",1,1," << rtol << ',' << StateName(result.state_) << ','
            << result.stats_.accepted_ << ',' << result.stats_.rejected_ << ',' << result.stats_.constraint_init_iterations_
            << ',' << RelativeError(state.variables_[0][x], exact) << ',' << maximum_algebraic_error << ','
            << maximum_residual << '\n';
      }
    }

    // All but one cell are exactly stationary. A global RMS controller dilutes the
    // one active cell's error by the number of zero-error cells.
    const std::array<std::size_t, 4> cell_counts{ 1, 10, 100, 1000 };
    for (const std::size_t cells : cell_counts)
    {
      auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      parameters.max_number_of_steps_ = 100000;
      auto solver = BuildSlavedSolver(parameters, 1);
      auto state = solver.GetState(cells);
      SetConditions(state);
      const double rtol = 1.0e-6;
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-9));
      const auto x = state.variable_map_.at("X");
      const auto z = state.variable_map_.at("Z0");
      state.variables_[cells - 1][x] = 1.0;
      state.variables_[cells - 1][z] = 1.0;

      const auto result = solver.Solve(1.0, state, parameters);
      const double exact = std::exp(-1.0);
      const double x_error = RelativeError(state.variables_[cells - 1][x], exact);
      const double z_error = RelativeError(state.variables_[cells - 1][z], exact);
      const double residual = std::abs(state.variables_[cells - 1][z] - state.variables_[cells - 1][x]);
      csv << "cell_rms_dilution,1," << cells << ",1," << rtol << ',' << StateName(result.state_) << ','
          << result.stats_.accepted_ << ',' << result.stats_.rejected_ << ',' << result.stats_.constraint_init_iterations_
          << ',' << x_error << ',' << z_error << ',' << residual << '\n';
    }

    const std::array<std::size_t, 4> solve_calls{ 1, 10, 100, 1000 };
    for (const std::size_t calls : solve_calls)
    {
      auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      parameters.max_number_of_steps_ = 100000;
      auto solver = BuildManufacturedSolver(parameters, 1.0, 1.0);
      auto state = solver.GetState(1);
      SetConditions(state);
      const double rtol = 1.0e-6;
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-9));
      const auto x = state.variable_map_.at("X");
      const auto z = state.variable_map_.at("Z");
      state.variables_[0][x] = 1.0;
      state.variables_[0][z] = 1.0;

      std::uint64_t steps = 0;
      std::uint64_t rejected = 0;
      std::uint64_t init_iterations = 0;
      double maximum_residual = 0.0;
      micm::SolverState final_state = micm::SolverState::NotYetCalled;
      for (std::size_t call = 0; call < calls; ++call)
      {
        const auto result = solver.Solve(1.0 / static_cast<double>(calls), state, parameters);
        final_state = result.state_;
        steps += result.stats_.accepted_;
        rejected += result.stats_.rejected_;
        init_iterations += result.stats_.constraint_init_iterations_;
        maximum_residual =
            std::max(maximum_residual, std::abs(state.variables_[0][z] - state.variables_[0][x] * state.variables_[0][x]));
        if (result.state_ != micm::SolverState::Converged)
        {
          break;
        }
      }
      csv << "output_cadence,1,1," << calls << ',' << rtol << ',' << StateName(final_state) << ',' << steps << ','
          << rejected << ',' << init_iterations << ',' << RelativeError(state.variables_[0][x], 0.5) << ','
          << RelativeError(state.variables_[0][z], 0.25) << ',' << maximum_residual << '\n';
    }

    // X crosses zero under X'=-1. The wrapper clamps differential X after the
    // Rosenbrock solve; the algebraic copy exposes whether that clamp is reprojected.
    {
      auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      auto solver = BuildSlavedSolver(parameters, 1, true);
      auto state = solver.GetState(1);
      SetConditions(state);
      state.SetRelativeTolerance(1.0e-6);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e-9));
      const auto x = state.variable_map_.at("X");
      const auto z = state.variable_map_.at("Z0");
      state.variables_[0][x] = 0.5;
      state.variables_[0][z] = 0.5;
      const auto result = solver.Solve(1.0, state, parameters);
      csv << "clamp_activation,1,1,1,1e-6," << StateName(result.state_) << ',' << result.stats_.accepted_ << ','
          << result.stats_.rejected_ << ',' << result.stats_.constraint_init_iterations_ << ',' << state.variables_[0][x]
          << ',' << state.variables_[0][z] << ',' << std::abs(state.variables_[0][z] - state.variables_[0][x]) << '\n';
    }

    summary << "- Controller/cadence: completed algebraic-count and public-Solve cadence sweeps.\n";
  }

  void RunVanDerPolEpsilon(const std::filesystem::path& output_directory, std::ostream& summary)
  {
    std::ofstream csv(output_directory / "vanderpol_epsilon.csv");
    csv << std::scientific << std::setprecision(17);
    csv << "epsilon,state,accepted,rejected,x_final,z_final,x_diff_vs_dae,z_diff_vs_dae\n";

    // Lienard-plane Van der Pol (negated fast variable, see the model classes)
    // from X0 = 2 with Z0 on the eps = 0 manifold (Z = X/(X^2-1)), integrated
    // to T = 0.5 — safely before the fold of the reduced problem at
    // t* = 3/2 - ln 2 ~ 0.807 where Z blows up at X = 1.
    constexpr double X0 = 2.0;
    constexpr double Z0 = X0 / (X0 * X0 - 1.0);
    constexpr double T_FINAL = 0.5;
    const double rtol = 1.0e-10;
    const double atol = 1.0e-13;

    auto parameters = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    parameters.max_number_of_steps_ = 1000000;

    struct VdpResult
    {
      micm::SolverState state;
      std::uint64_t accepted;
      std::uint64_t rejected;
      double x;
      double z;
    };
    auto run = [&](auto& solver) -> VdpResult
    {
      auto state = solver.GetState(1);
      SetConditions(state);
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, atol));
      state.variables_[0][state.variable_map_.at("X")] = X0;
      state.variables_[0][state.variable_map_.at("Z")] = Z0;
      const auto result = solver.Solve(T_FINAL, state, parameters);
      return { result.state_,
               result.stats_.accepted_,
               result.stats_.rejected_,
               state.variables_[0][state.variable_map_.at("X")],
               state.variables_[0][state.variable_map_.at("Z")] };
    };

    auto dae_solver = BuildVanDerPolDaeSolver(parameters);
    const auto dae = run(dae_solver);
    csv << 0.0 << ',' << StateName(dae.state) << ',' << dae.accepted << ',' << dae.rejected << ',' << dae.x << ','
        << dae.z << ",0,0\n";

    std::vector<double> epsilons;
    std::vector<double> x_diffs;
    bool all_converged = dae.state == micm::SolverState::Converged;
    for (int exponent = 0; exponent >= -6; --exponent)
    {
      const double epsilon = std::pow(10.0, exponent);
      auto solver = BuildVanDerPolOdeSolver(parameters, epsilon);
      const auto ode = run(solver);
      const double x_diff = std::abs(ode.x - dae.x);
      const double z_diff = std::abs(ode.z - dae.z);
      csv << epsilon << ',' << StateName(ode.state) << ',' << ode.accepted << ',' << ode.rejected << ',' << ode.x << ','
          << ode.z << ',' << x_diff << ',' << z_diff << '\n';
      epsilons.push_back(epsilon);
      x_diffs.push_back(x_diff);
      all_converged = all_converged && ode.state == micm::SolverState::Converged;
    }

    // Observed order of |x_ODE(eps) - x_DAE| over the two smallest decades;
    // singular-perturbation theory predicts first order in eps on the stable
    // slow manifold.
    const std::size_t n = x_diffs.size();
    const double observed_order =
        std::log10(x_diffs[n - 2] / x_diffs[n - 1]) / std::log10(epsilons[n - 2] / epsilons[n - 1]);
    summary << "- Van der Pol eps-limit: all runs " << (all_converged ? "converged" : "DID NOT ALL CONVERGE") << "; |x_ODE(eps)-x_DAE| = " << x_diffs.front()
            << " at eps=1 down to " << x_diffs.back() << " at eps=1e-6; observed order "
            << observed_order << " (expect ~1); DAE accepted steps " << dae.accepted << ".\n";
  }
}  // namespace

int main(int argc, char** argv)
{
  const std::filesystem::path output_directory = argc > 1 ? argv[1] : "dae_experiments";
  std::filesystem::create_directories(output_directory);

  std::ofstream summary_file(output_directory / "summary.md");
  summary_file << "# DAE convergence experiment run\n\n";
  summary_file << std::scientific << std::setprecision(6);

  std::cout << "Writing DAE convergence experiments to " << output_directory << '\n';
  RunProjectionScaling(output_directory, summary_file);
  RunFixedStepOrder(output_directory, summary_file);
  RunToleranceCoupling(output_directory, summary_file);
  RunNonlinearInitialization(output_directory, summary_file);
  RunControllerAndCadence(output_directory, summary_file);
  RunVanDerPolEpsilon(output_directory, summary_file);
  std::cout << "Completed DAE convergence experiments\n";
  return 0;
}
