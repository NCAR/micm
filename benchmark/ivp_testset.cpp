// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Community IVP test-set problems (CWI/Bari) run through MICM and judged
// against external Radau references that reproduce the published terminal
// values:
//
//   Chemical Akzo Nobel — index-1 DAE of dimension 6, M = diag(1,1,1,1,1,0),
//     0 = Ks*y1*y4 - y6. The kinetics have fractional reaction orders
//     (y1^4*sqrt(y2)), so the problem enters MICM as an external model.
//     Run twice: as the DAE, and as the exact reduced 5-species ODE with
//     y6 = Ks*y1*y4 substituted (the two must agree with each other and the
//     reference).
//
//   HIRES — stiff ODE of dimension 8 (plant physiology), t_end = 321.8122;
//     validates the base integrator on a published stiff problem.
#include <micm/CPU.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef MICM_BENCHMARK_DATA_DIR
#define MICM_BENCHMARK_DATA_DIR "benchmark/data"
#endif

namespace
{
  // Akzo Nobel constants (IVP test-set description).
  constexpr double AK_K1 = 18.7, AK_K2 = 0.58, AK_K3 = 0.09, AK_K4 = 0.42;
  constexpr double AK_BIGK = 34.4, AK_KLA = 3.3, AK_KS = 115.83;
  constexpr double AK_PCO2 = 0.9, AK_H = 737.0;

  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    const double lo = std::log10(t_first), hi = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, lo + (hi - lo) * i / (n - 1)));
    return times;
  }

  /// Chemical Akzo Nobel as the index-1 DAE: Y1..Y5 differential, Y6 algebraic.
  class AkzoDaeModel
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "Y1", "Y2", "Y3", "Y4", "Y5", "Y6" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5"), y6 = s.at("Y6");
      return { { y1, y1 }, { y1, y2 }, { y1, y3 }, { y1, y4 }, { y1, y5 },
               { y2, y1 }, { y2, y2 }, { y2, y4 }, { y2, y6 },
               { y3, y1 }, { y3, y2 }, { y3, y3 }, { y3, y4 }, { y3, y5 },
               { y4, y1 }, { y4, y3 }, { y4, y4 }, { y4, y5 },
               { y5, y1 }, { y5, y2 }, { y5, y3 }, { y5, y4 }, { y5, y5 }, { y5, y6 } };
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5"), y6 = s.at("Y6");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double v1 = state[cell][y1], v2 = state[cell][y2], v3 = state[cell][y3];
          const double v4 = state[cell][y4], v5 = state[cell][y5], v6 = state[cell][y6];
          const double sq2 = std::sqrt(std::max(v2, 0.0));
          const double r1 = AK_K1 * v1 * v1 * v1 * v1 * sq2;
          const double r2 = AK_K2 * v3 * v4;
          const double r3 = (AK_K2 / AK_BIGK) * v1 * v5;
          const double r4 = AK_K3 * v1 * v4 * v4;
          const double r5 = AK_K4 * v6 * v6 * sq2;
          const double f_in = AK_KLA * (AK_PCO2 / AK_H - v2);
          forcing[cell][y1] += -2.0 * r1 + r2 - r3 - r4;
          forcing[cell][y2] += -0.5 * r1 - r4 - 0.5 * r5 + f_in;
          forcing[cell][y3] += r1 - r2 + r3;
          forcing[cell][y4] += -r2 + r3 - 2.0 * r4;
          forcing[cell][y5] += r2 - r3 + r5;
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5"), y6 = s.at("Y6");
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double v1 = state[block][y1], v2 = state[block][y2], v3 = state[block][y3];
          const double v4 = state[block][y4], v5 = state[block][y5], v6 = state[block][y6];
          const double sq2 = std::sqrt(std::max(v2, 1.0e-300));
          const double dr1_1 = 4.0 * AK_K1 * v1 * v1 * v1 * sq2;
          const double dr1_2 = AK_K1 * v1 * v1 * v1 * v1 / (2.0 * sq2);
          const double dr2_3 = AK_K2 * v4;
          const double dr2_4 = AK_K2 * v3;
          const double dr3_1 = (AK_K2 / AK_BIGK) * v5;
          const double dr3_5 = (AK_K2 / AK_BIGK) * v1;
          const double dr4_1 = AK_K3 * v4 * v4;
          const double dr4_4 = 2.0 * AK_K3 * v1 * v4;
          const double dr5_2 = AK_K4 * v6 * v6 / (2.0 * sq2);
          const double dr5_6 = 2.0 * AK_K4 * v6 * sq2;
          jacobian[block][y1][y1] -= -2.0 * dr1_1 - dr3_1 - dr4_1;
          jacobian[block][y1][y2] -= -2.0 * dr1_2;
          jacobian[block][y1][y3] -= dr2_3;
          jacobian[block][y1][y4] -= dr2_4 - dr4_4;
          jacobian[block][y1][y5] -= -dr3_5;
          jacobian[block][y2][y1] -= -0.5 * dr1_1 - dr4_1;
          jacobian[block][y2][y2] -= -0.5 * dr1_2 - 0.5 * dr5_2 - AK_KLA;
          jacobian[block][y2][y4] -= -dr4_4;
          jacobian[block][y2][y6] -= -0.5 * dr5_6;
          jacobian[block][y3][y1] -= dr1_1 + dr3_1;
          jacobian[block][y3][y2] -= dr1_2;
          jacobian[block][y3][y3] -= -dr2_3;
          jacobian[block][y3][y4] -= -dr2_4;
          jacobian[block][y3][y5] -= dr3_5;
          jacobian[block][y4][y1] -= dr3_1 - 2.0 * dr4_1;
          jacobian[block][y4][y3] -= -dr2_3;
          jacobian[block][y4][y4] -= -dr2_4 - 2.0 * dr4_4;
          jacobian[block][y4][y5] -= dr3_5;
          jacobian[block][y5][y1] -= -dr3_1;
          jacobian[block][y5][y2] -= dr5_2;
          jacobian[block][y5][y3] -= dr2_3;
          jacobian[block][y5][y4] -= dr2_4;
          jacobian[block][y5][y5] -= -dr3_5;
          jacobian[block][y5][y6] -= dr5_6;
        }
      };
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "Y6" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "Y1", "Y4", "Y6" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y6 = s.at("Y6");
      return { { y6, s.at("Y1") }, { y6, s.at("Y4") }, { y6, y6 } };
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y1 = s.at("Y1"), y4 = s.at("Y4"), y6 = s.at("Y6");
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          forcing[cell][y6] = AK_KS * state[cell][y1] * state[cell][y4] - state[cell][y6];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto y1 = s.at("Y1"), y4 = s.at("Y4"), y6 = s.at("Y6");
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          jacobian[block][y6][y1] -= AK_KS * state[block][y4];
          jacobian[block][y6][y4] -= AK_KS * state[block][y1];
          jacobian[block][y6][y6] -= -1.0;
        }
      };
    }
  };

  /// The exact reduced 5-species ODE with y6 = Ks*y1*y4 substituted into r5.
  class AkzoReducedOdeModel
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "Y1", "Y2", "Y3", "Y4", "Y5" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5");
      std::set<std::pair<std::size_t, std::size_t>> elements;
      for (const auto row : { y1, y2, y3, y4, y5 })
        for (const auto col : { y1, y2, y3, y4, y5 })
          elements.insert({ row, col });
      return elements;
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double v1 = state[cell][y1], v2 = state[cell][y2], v3 = state[cell][y3];
          const double v4 = state[cell][y4], v5 = state[cell][y5];
          const double v6 = AK_KS * v1 * v4;
          const double sq2 = std::sqrt(std::max(v2, 0.0));
          const double r1 = AK_K1 * v1 * v1 * v1 * v1 * sq2;
          const double r2 = AK_K2 * v3 * v4;
          const double r3 = (AK_K2 / AK_BIGK) * v1 * v5;
          const double r4 = AK_K3 * v1 * v4 * v4;
          const double r5 = AK_K4 * v6 * v6 * sq2;
          const double f_in = AK_KLA * (AK_PCO2 / AK_H - v2);
          forcing[cell][y1] += -2.0 * r1 + r2 - r3 - r4;
          forcing[cell][y2] += -0.5 * r1 - r4 - 0.5 * r5 + f_in;
          forcing[cell][y3] += r1 - r2 + r3;
          forcing[cell][y4] += -r2 + r3 - 2.0 * r4;
          forcing[cell][y5] += r2 - r3 + r5;
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto y1 = s.at("Y1"), y2 = s.at("Y2"), y3 = s.at("Y3"), y4 = s.at("Y4"), y5 = s.at("Y5");
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double v1 = state[block][y1], v2 = state[block][y2], v3 = state[block][y3];
          const double v4 = state[block][y4], v5 = state[block][y5];
          const double v6 = AK_KS * v1 * v4;
          const double sq2 = std::sqrt(std::max(v2, 1.0e-300));
          const double dr1_1 = 4.0 * AK_K1 * v1 * v1 * v1 * sq2;
          const double dr1_2 = AK_K1 * v1 * v1 * v1 * v1 / (2.0 * sq2);
          const double dr2_3 = AK_K2 * v4;
          const double dr2_4 = AK_K2 * v3;
          const double dr3_1 = (AK_K2 / AK_BIGK) * v5;
          const double dr3_5 = (AK_K2 / AK_BIGK) * v1;
          const double dr4_1 = AK_K3 * v4 * v4;
          const double dr4_4 = 2.0 * AK_K3 * v1 * v4;
          // r5 = k4 * Ks^2 * y1^2 * y4^2 * sqrt(y2)
          const double dr5_1 = 2.0 * AK_K4 * AK_KS * AK_KS * v1 * v4 * v4 * sq2;
          const double dr5_2 = AK_K4 * v6 * v6 / (2.0 * sq2);
          const double dr5_4 = 2.0 * AK_K4 * AK_KS * AK_KS * v1 * v1 * v4 * sq2;
          jacobian[block][y1][y1] -= -2.0 * dr1_1 - dr3_1 - dr4_1;
          jacobian[block][y1][y2] -= -2.0 * dr1_2;
          jacobian[block][y1][y3] -= dr2_3;
          jacobian[block][y1][y4] -= dr2_4 - dr4_4;
          jacobian[block][y1][y5] -= -dr3_5;
          jacobian[block][y2][y1] -= -0.5 * dr1_1 - dr4_1 - 0.5 * dr5_1;
          jacobian[block][y2][y2] -= -0.5 * dr1_2 - 0.5 * dr5_2 - AK_KLA;
          jacobian[block][y2][y4] -= -dr4_4 - 0.5 * dr5_4;
          jacobian[block][y3][y1] -= dr1_1 + dr3_1;
          jacobian[block][y3][y2] -= dr1_2;
          jacobian[block][y3][y3] -= -dr2_3;
          jacobian[block][y3][y4] -= -dr2_4;
          jacobian[block][y3][y5] -= dr3_5;
          jacobian[block][y4][y1] -= dr3_1 - 2.0 * dr4_1;
          jacobian[block][y4][y3] -= -dr2_3;
          jacobian[block][y4][y4] -= -dr2_4 - 2.0 * dr4_4;
          jacobian[block][y4][y5] -= dr3_5;
          jacobian[block][y5][y1] -= -dr3_1 + dr5_1;
          jacobian[block][y5][y2] -= dr5_2;
          jacobian[block][y5][y3] -= dr2_3;
          jacobian[block][y5][y4] -= dr2_4 + dr5_4;
          jacobian[block][y5][y5] -= -dr3_5;
        }
      };
    }
  };

  /// HIRES stiff ODE (8 species).
  class HiresModel
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto h1 = s.at("H1"), h2 = s.at("H2"), h3 = s.at("H3"), h4 = s.at("H4");
      const auto h5 = s.at("H5"), h6 = s.at("H6"), h7 = s.at("H7"), h8 = s.at("H8");
      return { { h1, h1 }, { h1, h2 }, { h1, h3 },
               { h2, h1 }, { h2, h2 },
               { h3, h3 }, { h3, h4 }, { h3, h5 },
               { h4, h2 }, { h4, h3 }, { h4, h4 },
               { h5, h5 }, { h5, h6 }, { h5, h7 },
               { h6, h4 }, { h6, h5 }, { h6, h6 }, { h6, h7 }, { h6, h8 },
               { h7, h6 }, { h7, h7 }, { h7, h8 },
               { h8, h6 }, { h8, h7 }, { h8, h8 } };
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
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto h1 = s.at("H1"), h2 = s.at("H2"), h3 = s.at("H3"), h4 = s.at("H4");
      const auto h5 = s.at("H5"), h6 = s.at("H6"), h7 = s.at("H7"), h8 = s.at("H8");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double v1 = state[cell][h1], v2 = state[cell][h2], v3 = state[cell][h3], v4 = state[cell][h4];
          const double v5 = state[cell][h5], v6 = state[cell][h6], v7 = state[cell][h7], v8 = state[cell][h8];
          forcing[cell][h1] += -1.71 * v1 + 0.43 * v2 + 8.32 * v3 + 0.0007;
          forcing[cell][h2] += 1.71 * v1 - 8.75 * v2;
          forcing[cell][h3] += -10.03 * v3 + 0.43 * v4 + 0.035 * v5;
          forcing[cell][h4] += 8.32 * v2 + 1.71 * v3 - 1.12 * v4;
          forcing[cell][h5] += -1.745 * v5 + 0.43 * v6 + 0.43 * v7;
          forcing[cell][h6] += -280.0 * v6 * v8 + 0.69 * v4 + 1.71 * v5 - 0.43 * v6 + 0.69 * v7;
          forcing[cell][h7] += 280.0 * v6 * v8 - 1.81 * v7;
          forcing[cell][h8] += -280.0 * v6 * v8 + 1.81 * v7;
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto h1 = s.at("H1"), h2 = s.at("H2"), h3 = s.at("H3"), h4 = s.at("H4");
      const auto h5 = s.at("H5"), h6 = s.at("H6"), h7 = s.at("H7"), h8 = s.at("H8");
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double v6 = state[block][h6], v8 = state[block][h8];
          jacobian[block][h1][h1] -= -1.71;
          jacobian[block][h1][h2] -= 0.43;
          jacobian[block][h1][h3] -= 8.32;
          jacobian[block][h2][h1] -= 1.71;
          jacobian[block][h2][h2] -= -8.75;
          jacobian[block][h3][h3] -= -10.03;
          jacobian[block][h3][h4] -= 0.43;
          jacobian[block][h3][h5] -= 0.035;
          jacobian[block][h4][h2] -= 8.32;
          jacobian[block][h4][h3] -= 1.71;
          jacobian[block][h4][h4] -= -1.12;
          jacobian[block][h5][h5] -= -1.745;
          jacobian[block][h5][h6] -= 0.43;
          jacobian[block][h5][h7] -= 0.43;
          jacobian[block][h6][h4] -= 0.69;
          jacobian[block][h6][h5] -= 1.71;
          jacobian[block][h6][h6] -= -280.0 * v8 - 0.43;
          jacobian[block][h6][h7] -= 0.69;
          jacobian[block][h6][h8] -= -280.0 * v6;
          jacobian[block][h7][h6] -= 280.0 * v8;
          jacobian[block][h7][h7] -= -1.81;
          jacobian[block][h7][h8] -= 280.0 * v6;
          jacobian[block][h8][h6] -= -280.0 * v8;
          jacobian[block][h8][h7] -= 1.81;
          jacobian[block][h8][h8] -= -280.0 * v6;
        }
      };
    }
  };

  // Reference loading: header time,<name1>,<name2>,... -> map by column name.
  struct Reference
  {
    std::vector<double> times;
    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;  // one vector per name
  };

  Reference LoadReference(const std::string& path)
  {
    Reference ref;
    std::ifstream file(path);
    if (!file)
      return ref;
    std::string line;
    std::getline(file, line);
    {
      std::istringstream header(line);
      std::string field;
      std::getline(header, field, ',');  // time
      while (std::getline(header, field, ','))
        ref.names.push_back(field);
    }
    ref.columns.resize(ref.names.size());
    while (std::getline(file, line))
    {
      std::istringstream row(line);
      std::string field;
      std::vector<double> values;
      while (std::getline(row, field, ','))
        values.push_back(std::stod(field));
      if (values.size() != ref.names.size() + 1)
        continue;
      ref.times.push_back(values[0]);
      for (std::size_t i = 0; i < ref.names.size(); ++i)
        ref.columns[i].push_back(values[i + 1]);
    }
    return ref;
  }

  struct RunOutcome
  {
    std::uint64_t accepted = 0;
    std::uint64_t rejected = 0;
    bool converged = true;
    double max_rel_err = 0.0;
  };

  // Integrate `solver` over the output grid, comparing `species` (subset of the
  // reference columns, matched by name) at every output time.
  template<typename Solver>
  RunOutcome RunAgainstReference(
      Solver& solver,
      const std::vector<std::pair<std::string, double>>& initial_values,
      double rtol,
      double atol,
      const std::vector<double>& output_times,
      const Reference& ref,
      const std::vector<std::string>& species)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, atol));
    for (const auto& [name, value] : initial_values)
      state.variables_[0][state.variable_map_.at(name)] = value;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    std::vector<std::size_t> ref_columns;
    for (const auto& name : species)
    {
      for (std::size_t i = 0; i < ref.names.size(); ++i)
        if (ref.names[i] == name)
          ref_columns.push_back(i);
    }

    RunOutcome out;
    auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
    double current = 0.0;
    for (std::size_t i_out = 0; i_out < output_times.size(); ++i_out)
    {
      const double t_out = output_times[i_out];
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
        {
          out.converged = false;
          return out;
        }
        if (result.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          return out;
        }
        out.accepted += result.stats_.accepted_;
        out.rejected += result.stats_.rejected_;
        done += result.stats_.final_time_;
      }
      current = t_out;
      for (std::size_t k = 0; k < species.size(); ++k)
      {
        const double got = state.variables_[0][state.variable_map_.at(species[k])];
        out.max_rel_err = std::max(out.max_rel_err, rel(got, ref.columns[ref_columns[k]][i_out]));
      }
    }
    return out;
  }
}  // namespace

int main()
{
  const std::vector<double> rtols = { 1e-4, 1e-6, 1e-8, 1e-10 };
  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

  const auto akzo_ref = LoadReference(std::string(MICM_BENCHMARK_DATA_DIR) + "/akzo_reference.csv");
  const auto hires_ref = LoadReference(std::string(MICM_BENCHMARK_DATA_DIR) + "/hires_reference.csv");
  if (akzo_ref.times.empty() || hires_ref.times.empty())
  {
    std::cout << "missing references in " << MICM_BENCHMARK_DATA_DIR
              << "; generate with benchmark/generate_reference_solutions.py\n";
    return 1;
  }

  std::ofstream csv("ivp_testset.csv");
  csv << "problem,method,rtol,accepted,rejected,converged,max_rel_err\n";
  auto write_row = [&](const std::string& problem, const std::string& method, double rtol, const RunOutcome& r)
  {
    csv << problem << ',' << method << ',' << rtol << ',' << r.accepted << ',' << r.rejected << ','
        << (r.converged ? 1 : 0) << ',' << r.max_rel_err << '\n';
    std::cout << problem << " " << method << " rtol=" << rtol << "  accepted=" << r.accepted
              << "  max_rel_err=" << r.max_rel_err << (r.converged ? "" : "  [FAILED]") << "\n";
  };

  // ---- Chemical Akzo Nobel: DAE and exact reduced ODE ----
  {
    const auto output_times = OutputTimes(1e-1, 180.0, 21);
    const double y6_0 = AK_KS * 0.444 * 0.007;
    const std::vector<std::pair<std::string, double>> dae_init = { { "Y1", 0.444 }, { "Y2", 0.00123 }, { "Y3", 0.0 },
                                                                   { "Y4", 0.007 }, { "Y5", 0.0 },     { "Y6", y6_0 } };
    const std::vector<std::pair<std::string, double>> ode_init = { { "Y1", 0.444 }, { "Y2", 0.00123 }, { "Y3", 0.0 },
                                                                   { "Y4", 0.007 }, { "Y5", 0.0 } };
    const std::vector<std::string> dae_species = { "Y1", "Y2", "Y3", "Y4", "Y5", "Y6" };
    const std::vector<std::string> ode_species = { "Y1", "Y2", "Y3", "Y4", "Y5" };

    for (double rtol : rtols)
    {
      {
        std::vector<micm::PhaseSpecies> species;
        for (const auto& name : dae_species)
          species.push_back(micm::Species(name));
        const micm::Phase gas{ "gas", species };
        auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                          .SetSystem(micm::System(gas))
                          .AddExternalModel(AkzoDaeModel())
                          .SetReorderState(false)
                          .Build();
        write_row("akzo", "dae", rtol, RunAgainstReference(solver, dae_init, rtol, 1e-14, output_times, akzo_ref, dae_species));
      }
      {
        std::vector<micm::PhaseSpecies> species;
        for (const auto& name : ode_species)
          species.push_back(micm::Species(name));
        const micm::Phase gas{ "gas", species };
        auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                          .SetSystem(micm::System(gas))
                          .AddExternalModel(AkzoReducedOdeModel())
                          .SetReorderState(false)
                          .Build();
        write_row(
            "akzo", "reduced_ode", rtol, RunAgainstReference(solver, ode_init, rtol, 1e-14, output_times, akzo_ref, ode_species));
      }
    }
  }

  // ---- HIRES ----
  {
    const auto output_times = OutputTimes(1e-1, 321.8122, 21);
    const std::vector<std::pair<std::string, double>> init = { { "H1", 1.0 }, { "H2", 0.0 }, { "H3", 0.0 },
                                                               { "H4", 0.0 }, { "H5", 0.0 }, { "H6", 0.0 },
                                                               { "H7", 0.0 }, { "H8", 0.0057 } };
    const std::vector<std::string> species = { "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8" };
    for (double rtol : rtols)
    {
      std::vector<micm::PhaseSpecies> phase_species;
      for (const auto& name : species)
        phase_species.push_back(micm::Species(name));
      const micm::Phase gas{ "gas", phase_species };
      auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .AddExternalModel(HiresModel())
                        .SetReorderState(false)
                        .Build();
      write_row("hires", "ode", rtol, RunAgainstReference(solver, init, rtol, 1e-14, output_times, hires_ref, species));
    }
  }

  std::cout << "wrote ivp_testset.csv\n";
  return 0;
}
